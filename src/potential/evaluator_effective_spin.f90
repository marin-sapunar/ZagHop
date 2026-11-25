!--------------------------------------------------------------------------------------------------
! MODULE: evaluator_effective_spin_mod
!> @brief Effective spin states built from a potential that gives every spin component.
!> @details
!! The effective_spin_evaluator holds the original evaluator (the source) and presents a model
!! with one state for each spin-free state. The spin-free Hamiltonian is that of the source. The
!! spin-orbit couplings of the model are built from those of the source (soc_model):
!! - projection - The Hamiltonian of the source restricted to the singlets and the triplet states
!!   T_k u, where u is a fixed combination of the triplet components: H_eff = P^T H P.
!! - rss - Each singlet-triplet coupling is the root-sum-square over the components, s_ik |h_ik|.
!!   Couplings between triplets are dropped.
!! - rms - As rss, divided by sqrt(3).
!! - sq3 - As rss, divided by 3 (the SO-SQ3 scaling of effective operators in the literature).
!!
!! The signs s_ik of rss, rms and sq3 follow the couplings by continuity (Granucci, Persico and
!! Spighi, J. Chem. Phys. 137, 22A501 (2012)). At the first point, s_ik is the sign of u.h_ik.
!! Afterwards, s_ik changes when h_ik changes direction by more than 90 degrees from the previous
!! point (passed to eval), so the coupling stays smooth when it passes through zero. Without
!! signed, all s_ik are +1.
!!
!! u is the eigenvector with the largest eigenvalue of G = sum_ik h_ik h_ik^T, where h_ik are the
!! couplings of singlet i with the three components of triplet k, at a reference geometry. For
!! couplings of rank 1 (all h_ik parallel), projection and signed rss are exact, and rms and sq3
!! differ from them only by their scaling. The retained fraction, sum_ik (u.h_ik)**2 /
!! sum_ik |h_ik|**2, measures how far the couplings are from rank 1.
!!
!! Geometries, coordinate conversion, frames and atoms are delegated to the source. The source
!! only evaluates its Hamiltonian (eval_hamiltonian). Its bases are never built.
!!
!! The model is always built from all states of the source. If the dynamics of the source runs
!! on a subset of its states, the dynamics of the model runs on the same subset.
!--------------------------------------------------------------------------------------------------
module evaluator_effective_spin_mod
    use global_defs
    use evaluator_base_mod
    use evaluator_diabatic_mod
    implicit none

    private
    public :: effective_spin_evaluator


    real(dp), parameter :: ref_coupling_tol = 1.0e-10_dp !< Smallest root-sum-square of all
                                                         !! singlet-triplet couplings at the
                                                         !! reference geometry (Eh).
    real(dp), parameter :: degenerate_tol = 0.01_dp !< Relative difference of the two largest
                                                    !! eigenvalues of G below which u is not well
                                                    !! defined.


    type, extends(diabatic_evaluator) :: effective_spin_evaluator
        class(diabatic_evaluator), allocatable :: source !< Potential with all spin components.
        character(len=:), allocatable :: soc_model !< 'projection', 'rss', 'rms' or 'sq3'.
        logical :: signed = .true. !< Signs of rss, rms and sq3 by continuity, or all +1.
        real(dp) :: u(3) = 0.0_dp !< Combination of the triplet components.
        real(dp), allocatable :: proj(:, :) !< Isometry P from the states of the model to the
                                            !! states of the source, for all states of both.
        integer :: ms = 0 !< Manifold of the singlets (0 without singlets).
        integer :: mt = 0 !< Manifold of the triplets.
        real(dp) :: ref_fraction = 1.0_dp !< Retained fraction at the reference geometry.
        real(dp) :: fraction = 1.0_dp !< Retained fraction at the current geometry.
        real(dp) :: min_fraction = 1.0_dp !< Smallest retained fraction of all evaluations.
        real(dp), allocatable :: h_ref(:, :, :) !< Reference singlet-triplet couplings h(c, i, k)
        !! of the source for the signs at the next point (rss, rms, sq3): those of this point, or,
        !! for a pair whose couplings are zero here, its last nonzero couplings at earlier points.
        !! During eval, the references of the previous point.
        real(dp), allocatable :: sign_ref(:, :) !< Signs s_ik that belong to h_ref.
    contains
        procedure :: init => effective_init
        procedure :: eval => effective_eval
        procedure :: update_geometry => effective_update_geometry
        procedure :: eval_diab => effective_eval_diab
        procedure :: get_derivative_target_coord => effective_get_derivative_target_coord
        procedure :: get_reference_geometry => effective_get_reference_geometry
        procedure :: set_cartesian_frame => effective_set_cartesian_frame
        procedure :: check_atoms => effective_check_atoms
        procedure :: print_summary => effective_print_summary
        procedure :: add_state_layout => effective_add_state_layout
        procedure :: diabatic_from_source => effective_diabatic_from_source
        procedure :: save_state => effective_save_state
        procedure :: load_state => effective_load_state
        procedure :: st_couplings
        procedure :: retained_fraction
        procedure :: set_rss_couplings
        procedure :: coupling_signs
    end type effective_spin_evaluator


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_init
    !> @brief Build the effective model from a potential with all spin components.
    !> @details
    !! A copy of the source is kept. Its spin manifolds must be set, with spin_states =
    !! components, and only singlets and triplets. u is computed at the reference geometry, given
    !! in ref_geom or, by default, the default_reference of the source. Nothing of the evaluation at
    !! the reference is kept: the first point of the trajectory overwrites the state of the
    !! source.
    !----------------------------------------------------------------------------------------------
    subroutine effective_init(self, src, soc_model, ref_geom, ref_coordinates, signed)
        use linalg_wrapper_mod, only : syev
        class(effective_spin_evaluator), intent(inout) :: self
        class(potential_evaluator), intent(in) :: src !< Potential with all spin components.
        character(len=*), intent(in) :: soc_model !< 'projection', 'rss', 'rms' or 'sq3'.
        real(dp), intent(in), optional :: ref_geom(:, :) !< Reference geometry.
        character(len=*), intent(in), optional :: ref_coordinates !< Coordinates of ref_geom.
        logical, intent(in), optional :: signed !< Signs of rss, rms and sq3 by continuity.
        real(dp), allocatable :: geom(:, :)
        real(dp), allocatable :: h(:, :, :)
        real(dp) :: g(3, 3), lambda(3), tr_g
        integer :: m, k, c

        select case(soc_model)
        case('projection', 'rss', 'rms', 'sq3')
            self%soc_model = soc_model
        case default
            call init_error('Unrecognized effective spin-orbit coupling model: '//soc_model//     &
            &               '. Expected projection, rss, rms or sq3.')
        end select
        self%signed = .true.
        if (present(signed)) self%signed = signed
        if (soc_model == 'projection' .and. .not. self%signed) then
            call init_error('The projection model is always signed. unsigned is only '//          &
            &               'available for rss, rms and sq3.')
        end if

        select type(src)
        class is (diabatic_evaluator)
            allocate(self%source, source=src)
        class default
            call init_error('Effective spin states need a potential given by a diabatic '//       &
            &               'Hamiltonian.')
        end select

        ! Spin manifolds of the source.
        if (.not. self%source%model_states%is_set()) then
            call init_error('The spin manifolds of the potential are not set.')
        end if
        if (self%source%model_states%spin_states /= 'components') then
            call init_error('Effective spin states need a potential with spin_states = '//        &
            &               'components.')
        end if
        if (any(self%source%model_states%manifold%mult /= 1 .and.                                &
        &       self%source%model_states%manifold%mult /= 3)) then
            call init_error('Effective spin states are only implemented for singlets and '//      &
            &               'triplets.')
        end if
        do m = 1, size(self%source%model_states%manifold)
            if (self%source%model_states%manifold(m)%mult == 1) self%ms = m
            if (self%source%model_states%manifold(m)%mult == 3) self%mt = m
        end do
        if (self%mt == 0) call init_error('Effective spin states need triplet states.')

        ! The model has the manifolds of the source, with one state for each spin-free state,
        ! and the dynamics runs on the same spin-free states as that of the source.
        self%classical_coordinates = self%source%classical_coordinates
        self%ndof = self%source%ndof
        call self%set_model_states(self%source%model_states%manifold%n,                           &
        &                          self%source%model_states%manifold%mult, 'effective_rss')
        call self%select_states(self%source%states%manifold%n, self%source%states%manifold%mult)
        self%exact_gradient = self%source%exact_gradient
        ! The source checks the couplings within its groups, in its own layout.
        self%soc_checked = .true.

        ! Couplings at the reference geometry.
        if (present(ref_geom) .neqv. present(ref_coordinates)) then
            call init_error('ref_geom and ref_coordinates are given together.')
        end if
        if (present(ref_geom)) then
            call self%source%update_geometry(ref_geom, ref_coordinates)
        else
            call self%source%default_reference(geom)
            if (.not. allocated(geom)) then
                call init_error('The potential has no default reference geometry. Give one '//    &
                &               'with reference_geometry.')
            end if
            call self%source%update_geometry(geom, self%source%classical_coordinates)
        end if
        call self%source%eval_hamiltonian()

        ! u is the dominant direction of the couplings in the space of the triplet components.
        h = self%st_couplings(self%source%soc_w)
        g = matmul(reshape(h, [3, size(h) / 3]), transpose(reshape(h, [3, size(h) / 3])))
        tr_g = g(1, 1) + g(2, 2) + g(3, 3)
        if (sqrt(tr_g) < ref_coupling_tol) then
            call init_error('No singlet-triplet spin-orbit coupling at the reference geometry. '// &
            &               'Give another one with reference_geometry.')
        end if
        call syev(g, lambda, jobz='V', uplo='U')
        self%u = g(:, 3)
        if (self%u(maxloc(abs(self%u), dim=1)) < 0.0_dp) self%u = -self%u
        self%ref_fraction = lambda(3) / tr_g
        if (lambda(3) - lambda(2) < degenerate_tol * lambda(3)) then
            write(stderr, *) 'Warning in evaluator_effective_spin_mod, effective_init subroutine.'
            write(stderr, *) '  The two largest eigenvalues of the singlet-triplet coupling ',     &
            &                'matrix G at the reference geometry'
            write(stderr, '(3x,a,2e12.4,a)') 'are nearly degenerate (', lambda(3:2:-1),           &
            &                                ' Eh^2). The combination u of the triplet'
            write(stderr, *) '  components is not well defined.'
        end if

        ! Isometry from the states of the model to the states of the source.
        allocate(self%proj(self%source%model_states%n_state, self%model_states%n_sf),             &
        &        source=0.0_dp)
        do m = 1, size(self%source%model_states%manifold)
            associate(mf => self%source%model_states%manifold(m))
                do c = 1, mf%n_comp
                    do k = 1, mf%n
                        if (mf%n_comp == 1) then
                            self%proj(mf%i0 + k - 1, mf%sf0 + k - 1) = 1.0_dp
                        else
                            self%proj(mf%i0 + (c - 1) * mf%n + k - 1, mf%sf0 + k - 1) = self%u(c)
                        end if
                    end do
                end do
            end associate
        end do

        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Effective spin states built from the spin components '//     &
            &                       '(soc_model = '//self%soc_model//').'
            if (self%soc_model /= 'projection') then
                if (self%signed) then
                    write(stdout, '(3x,a)') 'Signs of the singlet-triplet couplings by continuity.'
                else
                    write(stdout, '(3x,a)') 'Unsigned singlet-triplet couplings.'
                end if
            end if
            write(stdout, '(3x,a,3f10.6)') 'Combination of the triplet components, u:', self%u
            write(stdout, '(3x,a,f10.6)') 'Retained fraction of the singlet-triplet couplings '// &
            &                             'at the reference:', self%ref_fraction
        end if

    contains

        subroutine init_error(msg)
            character(len=*), intent(in) :: msg

            write(stderr, *) 'Error in evaluator_effective_spin_mod, effective_init subroutine.'
            write(stderr, *) '  '//msg
            stop
        end subroutine init_error

    end subroutine effective_init


    !> Pass the geometry to the source.
    subroutine effective_update_geometry(self, geometry, coordinates)
        class(effective_spin_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        call self%source%update_geometry(geometry, coordinates)
    end subroutine effective_update_geometry


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_eval
    !> @brief Evaluate the model at the current geometry, with the signs of the couplings by
    !! continuity from the previous point prev.
    !> @details
    !! The sign references of prev (h_ref, sign_ref) are copied before the evaluation, and
    !! set_rss_couplings replaces them with those of this point. Without prev, the signs are those
    !! of a first point. The rest is done by eval_all.
    !----------------------------------------------------------------------------------------------
    subroutine effective_eval(self, prev)
        class(effective_spin_evaluator), intent(inout) :: self
        class(potential_evaluator), intent(in), optional :: prev

        if (allocated(self%h_ref)) deallocate(self%h_ref)
        if (allocated(self%sign_ref)) deallocate(self%sign_ref)
        if (present(prev)) then
            select type(prev)
            class is (effective_spin_evaluator)
                if (allocated(prev%h_ref)) then
                    self%h_ref = prev%h_ref
                    self%sign_ref = prev%sign_ref
                end if
            class default
                write(stderr, *) 'Error in evaluator_effective_spin_mod, effective_eval subroutine.'
                write(stderr, *) '  The previous point is not given by effective spin states.'
                stop
            end select
        end if
        call eval_all(self, prev)
    end subroutine effective_eval


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_eval_diab
    !> @brief Evaluate the source and build the spin-free Hamiltonian and the effective spin-orbit
    !! couplings of the model.
    !----------------------------------------------------------------------------------------------
    subroutine effective_eval_diab(self)
        class(effective_spin_evaluator), intent(inout) :: self
        integer :: imode

        call self%source%eval_hamiltonian()

        ! The spin-free Hamiltonian is the same in both layouts.
        self%sf_w = self%source%sf_w
        self%sf_dw = self%source%sf_dw

        select case(self%soc_model)
        case('projection')
            ! P^T expand(sf_w) P = sf_w, because |u| = 1, so only the couplings are projected.
            self%soc_w = matmul(transpose(self%proj), matmul(self%source%soc_w, self%proj))
            do imode = 1, self%ndof
                self%soc_dw(imode, :, :) = matmul(transpose(self%proj),                           &
                &                                 matmul(self%source%soc_dw(imode, :, :), self%proj))
            end do
        case('rss', 'rms', 'sq3')
            call self%set_rss_couplings()
        end select

        self%fraction = self%retained_fraction(self%st_couplings(self%source%soc_w))
        self%min_fraction = min(self%min_fraction, self%fraction)
        if (stdp3) then
            write(stdout, '(5x,a,f10.6)') 'Retained fraction of the singlet-triplet couplings:',  &
            &                             self%fraction
        end if
    end subroutine effective_eval_diab


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_rss_couplings
    !> @brief Singlet-triplet couplings of the rss, rms and sq3 models.
    !> @details
    !! The coupling of singlet i and triplet k is s_ik |h_ik| (coupling_signs), divided by sqrt(3)
    !! for rms and by 3 for sq3. Its derivative is s_ik h_ik.dh_ik / |h_ik|, and zero where
    !! h_ik = 0. All other couplings are zero.
    !!
    !! The couplings and signs of the pairs with nonzero couplings then replace the references
    !! (h_ref, sign_ref) for the signs at the next point. A pair whose couplings are zero keeps its
    !! last nonzero couplings, so a sign change exactly at a point is still found at the next one.
    !----------------------------------------------------------------------------------------------
    subroutine set_rss_couplings(self)
        class(effective_spin_evaluator), intent(inout) :: self
        real(dp), allocatable :: h(:, :, :), dh(:, :, :)
        real(dp), allocatable :: s(:, :), r(:, :)
        integer :: i, k, imode, i0, k0

        h = self%st_couplings(self%source%soc_w)
        r = norm2(h, dim=1)
        s = self%coupling_signs(h)
        if (.not. allocated(self%h_ref)) then
            allocate(self%h_ref(3, size(h, 2), size(h, 3)), source=0.0_dp)
            allocate(self%sign_ref(size(h, 2), size(h, 3)), source=1.0_dp)
        end if
        do k = 1, size(h, 3)
            do i = 1, size(h, 2)
                if (r(i, k) > tiny(1.0_dp)) then
                    self%h_ref(:, i, k) = h(:, i, k)
                    self%sign_ref(i, k) = s(i, k)
                end if
            end do
        end do
        select case(self%soc_model)
        case('rms')
            s = s / sqrt(3.0_dp)
        case('sq3')
            s = s / 3.0_dp
        end select

        self%soc_w = 0.0_dp
        self%soc_dw = 0.0_dp
        if (self%ms == 0) return
        i0 = self%model_states%manifold(self%ms)%i0
        k0 = self%model_states%manifold(self%mt)%i0
        do k = 1, size(h, 3)
            do i = 1, size(h, 2)
                self%soc_w(i0 + i - 1, k0 + k - 1) = s(i, k) * r(i, k)
                self%soc_w(k0 + k - 1, i0 + i - 1) = s(i, k) * r(i, k)
            end do
        end do
        do imode = 1, self%ndof
            dh = self%st_couplings(self%source%soc_dw(imode, :, :))
            do k = 1, size(h, 3)
                do i = 1, size(h, 2)
                    if (r(i, k) <= tiny(1.0_dp)) cycle
                    self%soc_dw(imode, i0 + i - 1, k0 + k - 1) = s(i, k)                          &
                    &   * dot_product(h(:, i, k), dh(:, i, k)) / r(i, k)
                    self%soc_dw(imode, k0 + k - 1, i0 + i - 1)                                    &
                    &   = self%soc_dw(imode, i0 + i - 1, k0 + k - 1)
                end do
            end do
        end do
    end subroutine set_rss_couplings


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: coupling_signs
    !> @brief Signs s_ik of the rss, rms and sq3 couplings at the current geometry.
    !> @details
    !! Without signed, all signs are +1. Otherwise, a pair with nonzero couplings at an earlier
    !! point keeps the sign of its reference (sign_ref, of the previous point during eval), and
    !! changes it if h_ik(ref).h_ik < 0 (Granucci, Persico and Spighi 2012, Eq. 16). A pair
    !! without earlier couplings, at the first point or after only zero couplings, takes the sign
    !! of u.h_ik (+1 at zero), which gives the relative signs of the projection model. The signs
    !! depend only on the references of the previous point, so a second evaluation at the same
    !! geometry with the same previous point gives the same signs.
    !----------------------------------------------------------------------------------------------
    pure function coupling_signs(self, h) result(s)
        class(effective_spin_evaluator), intent(in) :: self
        real(dp), intent(in) :: h(:, :, :) !< Couplings h(c, i, k), from st_couplings.
        real(dp), allocatable :: s(:, :)
        integer :: i, k

        allocate(s(size(h, 2), size(h, 3)), source=1.0_dp)
        if (.not. self%signed) return
        do k = 1, size(h, 3)
            do i = 1, size(h, 2)
                if (allocated(self%h_ref)) then
                    if (norm2(self%h_ref(:, i, k)) > tiny(1.0_dp)) then
                        s(i, k) = self%sign_ref(i, k)
                        if (dot_product(self%h_ref(:, i, k), h(:, i, k)) < 0.0_dp) then
                            s(i, k) = -s(i, k)
                        end if
                        cycle
                    end if
                end if
                if (dot_product(self%u, h(:, i, k)) < 0.0_dp) s(i, k) = -1.0_dp
            end do
        end do
    end function coupling_signs


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: st_couplings
    !> @brief Couplings h(c, i, k) between singlet i and component c of triplet k, from a matrix
    !! over the states of the source.
    !----------------------------------------------------------------------------------------------
    pure function st_couplings(self, soc) result(h)
        class(effective_spin_evaluator), intent(in) :: self
        real(dp), intent(in) :: soc(:, :) !< Matrix over the states of the source.
        real(dp), allocatable :: h(:, :, :)
        integer :: ns, nt, s0, t0, k, c

        nt = self%source%model_states%manifold(self%mt)%n
        t0 = self%source%model_states%manifold(self%mt)%i0
        ns = 0
        s0 = 1
        if (self%ms > 0) then
            ns = self%source%model_states%manifold(self%ms)%n
            s0 = self%source%model_states%manifold(self%ms)%i0
        end if
        allocate(h(3, ns, nt))
        do k = 1, nt
            do c = 1, 3
                h(c, :, k) = soc(s0:s0+ns-1, t0 + (c - 1) * nt + k - 1)
            end do
        end do
    end function st_couplings


    !> Fraction of the squared singlet-triplet couplings along u, sum (u.h)**2 / sum |h|**2. One
    !! without couplings.
    pure function retained_fraction(self, h) result(f)
        class(effective_spin_evaluator), intent(in) :: self
        real(dp), intent(in) :: h(:, :, :) !< Couplings h(c, i, k), from st_couplings.
        real(dp) :: f
        real(dp) :: total

        total = sum(h**2)
        if (total > 0.0_dp) then
            f = sum(matmul(self%u, reshape(h, [3, size(h) / 3]))**2) / total
        else
            f = 1.0_dp
        end if
    end function retained_fraction


    !> Derivatives are converted by the source, which holds the last geometry.
    function effective_get_derivative_target_coord(self, d) result(dr)
        class(effective_spin_evaluator), intent(in) :: self
        real(dp), intent(in) :: d(:)
        real(dp), allocatable :: dr(:)

        dr = self%source%get_derivative_target_coord(d)
    end function effective_get_derivative_target_coord


    subroutine effective_get_reference_geometry(self, coordinates, geom)
        class(effective_spin_evaluator), intent(in) :: self
        character(len=*), intent(in) :: coordinates
        real(dp), allocatable, intent(out) :: geom(:, :)

        call self%source%get_reference_geometry(coordinates, geom)
    end subroutine effective_get_reference_geometry


    subroutine effective_set_cartesian_frame(self, frame)
        class(effective_spin_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: frame

        call self%source%set_cartesian_frame(frame)
    end subroutine effective_set_cartesian_frame


    subroutine effective_check_atoms(self, sym, mass)
        class(effective_spin_evaluator), intent(in) :: self
        character(len=*), intent(in) :: sym(:)
        real(dp), intent(in) :: mass(:)

        call self%source%check_atoms(sym, mass)
    end subroutine effective_check_atoms


    !> Smallest retained fraction of the singlet-triplet couplings along the trajectory.
    subroutine effective_print_summary(self, unit)
        class(effective_spin_evaluator), intent(in) :: self
        integer, intent(in) :: unit

        write(unit, '(1x,a,f10.6)') 'Smallest retained fraction of the singlet-triplet '//        &
        &                           'couplings:', self%min_fraction
    end subroutine effective_print_summary


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_add_state_layout
    !> @brief Add the layout of the states of the model, the effective model and the layout of
    !! the source to a JSON object.
    !> @details
    !! effective_spin_states holds the model (soc_model), whether the couplings are signed, u and
    !! the retained fraction at the reference geometry. source holds the layout of the source.
    !----------------------------------------------------------------------------------------------
    subroutine effective_add_state_layout(self, json, obj)
        use json_module, only : json_core, json_value
        class(effective_spin_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the layout is added.
        type(json_value), pointer :: eff, src

        call add_state_layout(self, json, obj)
        call json%create_object(eff, 'effective_spin_states')
        call json%add(eff, 'model', self%soc_model)
        call json%add(eff, 'signed', self%signed)
        call json%add(eff, 'u', self%u)
        call json%add(eff, 'reference_fraction', self%ref_fraction)
        call json%add(obj, eff)
        call json%create_object(src, 'source')
        call self%source%add_state_layout(json, src)
        call json%add(obj, src)
    end subroutine effective_add_state_layout


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: effective_diabatic_from_source
    !> @brief Diabatic state of the model given by the index of a state of the source.
    !> @details
    !! The index follows the component_order of the source. Each spin component of the source is
    !! mapped to its spin-free state, which is the state of the model. The result is 0 if the
    !! source has no state with this index.
    !----------------------------------------------------------------------------------------------
    pure function effective_diabatic_from_source(self, i_source) result(state)
        class(effective_spin_evaluator), intent(in) :: self
        integer, intent(in) :: i_source !< Index of the state in the source.
        integer :: state
        type(state_location) :: loc

        state = self%source%diabatic_from_source(i_source)
        if (state > 0) then
            loc = self%source%model_states%locate_state(state)
            state = loc%sf
        end if
    end function effective_diabatic_from_source

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_save_state
    !> @brief Add the state of the evaluated point to a JSON object, with the effective model in
    !! the object 'effective_spin_states'.
    !> @details
    !! The state of the source (its geometry and Hamiltonian) is saved in the object 'source'.
    !! u and the projection are set up from the input, so only the quantities of the point (the
    !! retained fractions and the sign references) are saved.
    !----------------------------------------------------------------------------------------------
    subroutine effective_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        class(effective_spin_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: eff, src

        call diabatic_save_state(self, json, obj)
        call json%create_object(eff, 'effective_spin_states')
        call json%add(eff, 'soc_model', self%soc_model)
        call json%add(eff, 'fraction', self%fraction)
        call json%add(eff, 'min_fraction', self%min_fraction)
        if (allocated(self%h_ref)) then
            call json_add_array(json, eff, 'h_ref', self%h_ref)
            call json_add_array(json, eff, 'sign_ref', self%sign_ref)
        end if
        call json%create_object(src, 'source')
        call self%source%save_state(json, src)
        call json%add(eff, src)
        call json%add(obj, eff)
    end subroutine effective_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: effective_load_state
    !> @brief Restore the state of an evaluated point saved by save_state. Stops if the model of
    !! the couplings differs.
    !----------------------------------------------------------------------------------------------
    subroutine effective_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child, json_check_value
        class(effective_spin_evaluator), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: eff

        call diabatic_load_state(self, json, obj)
        eff => json_get_child(json, obj, 'effective_spin_states')
        call json_check_value(json, eff, 'soc_model', self%soc_model)
        call json_get_required(json, eff, 'fraction', self%fraction)
        call json_get_required(json, eff, 'min_fraction', self%min_fraction)
        if (allocated(self%h_ref)) deallocate(self%h_ref)
        if (allocated(self%sign_ref)) deallocate(self%sign_ref)
        if (json%valid_path(eff, 'h_ref')) then
            call json_get_required(json, eff, 'h_ref', self%h_ref)
            call json_get_required(json, eff, 'sign_ref', self%sign_ref,                          &
            &                      [size(self%h_ref, 2), size(self%h_ref, 3)])
        end if
        call self%source%load_state(json, json_get_child(json, eff, 'source'))
    end subroutine effective_load_state

end module evaluator_effective_spin_mod
