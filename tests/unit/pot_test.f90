module pot_test_mod
    use global_defs, only : dp
    use evaluator_base_mod, only : potential_evaluator
    use evaluator_diabatic_mod, only : diabatic_evaluator
    implicit none

    !> Evaluator without a potential, used to test the spin manifold setup.
    type, extends(diabatic_evaluator) :: group_evaluator
    contains
        procedure :: set_geometry => group_set_geometry
        procedure :: eval_diab => group_eval_diab
    end type group_evaluator

    !> Evaluator whose derivatives are converted to coordinates with three times as many
    !! elements, (d, 2 d, 3 d).
    type, extends(group_evaluator) :: coord_evaluator
    contains
        procedure :: get_derivative_target_coord => coord_derivative
    end type coord_evaluator

    !> Evaluator with a spin-free Hamiltonian and spin-orbit couplings linear in the coordinates,
    !! w = w0 + sum_n x(n) w1(n, :, :). The geometry is x, with dimensions (1, ndof), and the
    !! default reference geometry is x = 0. With scaled, derivatives are converted to coordinates
    !! with three times as many elements, (d, 2 d, 3 d).
    type, extends(diabatic_evaluator) :: linear_evaluator
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: sf0(:, :), sf1(:, :, :)
        real(dp), allocatable :: soc0(:, :), soc1(:, :, :)
        logical :: scaled = .false.
    contains
        procedure :: set_geometry => linear_set_geometry
        procedure :: eval_diab => linear_eval_diab
        procedure :: default_reference => linear_default_reference
        procedure :: get_derivative_target_coord => linear_derivative
    end type linear_evaluator

    !> Potential stored as in the trajectory data, to test copies.
    type :: pot_holder
        class(potential_evaluator), allocatable :: pot
    end type pot_holder

contains

    subroutine linear_set_geometry(self, geometry, coordinates)
        class(linear_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        self%x = geometry(1, :)
    end subroutine linear_set_geometry

    subroutine linear_eval_diab(self)
        class(linear_evaluator), intent(inout) :: self
        integer :: n

        self%sf_w = self%sf0
        self%soc_w = self%soc0
        do n = 1, self%ndof
            self%sf_w = self%sf_w + self%x(n) * self%sf1(n, :, :)
            self%soc_w = self%soc_w + self%x(n) * self%soc1(n, :, :)
        end do
        self%sf_dw = self%sf1
        self%soc_dw = self%soc1
    end subroutine linear_eval_diab

    subroutine linear_default_reference(self, geom)
        class(linear_evaluator), intent(in) :: self
        real(dp), allocatable, intent(out) :: geom(:, :)

        allocate(geom(1, self%ndof), source=0.0_dp)
    end subroutine linear_default_reference

    function linear_derivative(self, d) result(dr)
        class(linear_evaluator), intent(in) :: self
        real(dp), intent(in) :: d(:)
        real(dp), allocatable :: dr(:)

        if (self%scaled) then
            dr = [d, 2.0_dp * d, 3.0_dp * d]
        else
            dr = d
        end if
    end function linear_derivative

    subroutine group_set_geometry(self, geometry, coordinates)
        class(group_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
    end subroutine group_set_geometry

    subroutine group_eval_diab(self)
        class(group_evaluator), intent(inout) :: self
    end subroutine group_eval_diab

    function coord_derivative(self, d) result(dr)
        class(coord_evaluator), intent(in) :: self
        real(dp), intent(in) :: d(:)
        real(dp), allocatable :: dr(:)

        dr = [d, 2.0_dp * d, 3.0_dp * d]
    end function coord_derivative

end module pot_test_mod


program pot_test
    implicit none
    integer :: ret
    character(len=32) :: arg

    ! With an argument, run a case which should stop with an error.
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        call effective_error(trim(arg))
        write(*, *) 'No error.'
        stop 1
    end if

    ret = group_test()
    if (ret /= 0) then
        write(*,*) "Error in group_test."
        stop ret
    end if

    ret = vc_group_test()
    if (ret /= 0) then
        write(*,*) "Error in vc_group_test."
        stop ret
    end if

    ret = spin_free_test()
    if (ret /= 0) then
        write(*,*) "Error in spin_free_test."
        stop ret
    end if

    ret = component_order_test()
    if (ret /= 0) then
        write(*,*) "Error in component_order_test."
        stop ret
    end if

    ret = derivative_test()
    if (ret /= 0) then
        write(*,*) "Error in derivative_test."
        stop ret
    end if

    ret = effective_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_test."
        stop ret
    end if

    ret = effective_rotation_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_rotation_test."
        stop ret
    end if

    ret = effective_gradient_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_gradient_test."
        stop ret
    end if

    ret = effective_sign_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_sign_test."
        stop ret
    end if

    ret = effective_copy_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_copy_test."
        stop ret
    end if

    ret = model_1d_test()
    if (ret /= 0) then
        write(*,*) "Error in model_1d_test."
        stop ret
    end if

    ret = state_layout_test()
    if (ret /= 0) then
        write(*,*) "Error in state_layout_test."
        stop ret
    end if

    contains


    !> A full Hamiltonian given with all components of each state consecutive (as in Quantics) is
    !! reordered to the groups of components and split into the spin-free Hamiltonian and the
    !! spin-orbit couplings.
    integer function component_order_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        type(group_evaluator) :: pot, pot_zh, pot_eff
        real(dp) :: w(7, 7), dw(1, 7, 7), soc(7, 7), dsoc(7, 7)
        integer :: i

        ! One singlet and two triplets. Source order: S1, T1(1), T1(2), T1(3), T2(1), T2(2),
        ! T2(3). ZagHop order: S1, T1(1), T2(1), T1(2), T2(2), T1(3), T2(3).
        ret = 1
        pot%component_order = 'by_state'
        pot%ndof = 1
        call pot%set_spin_manifolds([1, 2], [1, 3], 'components')
        if (any(pot%state_map /= [1, 2, 5, 3, 6, 4, 7])) return
        call pot_zh%set_spin_manifolds([1, 2], [1, 3], 'components')
        if (any(pot_zh%state_map /= [(i, i = 1, 7)])) return
        pot_eff%component_order = 'by_state'
        call pot_eff%set_spin_manifolds([1, 2], [1, 3], 'effective_rss')
        if (any(pot_eff%state_map /= [1, 2, 3])) return

        ! Couplings of S1 with T1(1) and T2(2), and a derivative of the coupling with T2(2).
        ret = 2
        w = 0.0_dp
        dw = 0.0_dp
        do i = 2, 4
            w(i, i) = 0.1_dp
        end do
        do i = 5, 7
            w(i, i) = 0.2_dp
        end do
        w(1, 2) = 1.0e-4_dp
        w(1, 6) = 2.0e-4_dp
        dw(1, 1, 6) = 3.0e-4_dp
        w(:, 1) = w(1, :)
        dw(1, :, 1) = dw(1, 1, :)
        call pot%set_full_hamiltonian(w, dw)
        if (any(pot%sf_w /= reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp,    &
        &                            0.0_dp, 0.2_dp], [3, 3]))) return
        soc = 0.0_dp
        soc(1, 2) = 1.0e-4_dp
        soc(1, 5) = 2.0e-4_dp
        soc(:, 1) = soc(1, :)
        if (any(pot%soc_w /= soc)) return
        dsoc = 0.0_dp
        dsoc(1, 5) = 3.0e-4_dp
        dsoc(5, 1) = 3.0e-4_dp
        if (any(pot%soc_dw(1, :, :) /= dsoc)) return
        ret = 0
    end function component_order_test


    !> Gradients and coupling vectors are converted by get_derivative_target_coord and returned
    !! flat. update_geometry saves the eigenvectors of the previous geometry, and the spin-orbit
    !! couplings within groups are checked at the first evaluation.
    integer function derivative_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : coord_evaluator
        type(coord_evaluator) :: pot
        real(dp), allocatable :: grad(:), nadv(:), d(:), trans(:, :), sf_trans(:, :)
        real(dp) :: geom(3, 2)

        call pot%set_spin_manifolds([2], [1], 'components')
        pot%ndof = 2
        call pot%allocate_hamiltonian()
        pot%sf_w = reshape([0.0_dp, 0.01_dp, 0.01_dp, 0.1_dp], [2, 2])
        pot%sf_dw(1, :, :) = reshape([1.0_dp, 0.5_dp, 0.5_dp, -1.0_dp], [2, 2])
        pot%sf_dw(2, :, :) = reshape([0.2_dp, -0.3_dp, -0.3_dp, 0.4_dp], [2, 2])
        ret = 1
        if (pot%soc_checked) return
        call pot%eval()
        if (.not. pot%soc_checked) return

        ! Gradients are converted and flat, before and after a geometry is set.
        ret = 2
        grad = pot%get_gradient('diabatic', 1)
        d = pot%diab_dw(:, 1, 1)
        if (size(grad) /= 6) return
        if (any(grad /= [d, 2.0_dp * d, 3.0_dp * d])) return

        ret = 3
        trans = pot%adiab_trans
        sf_trans = pot%sf_trans
        geom = 0.0_dp
        call pot%update_geometry(geom, 'model')
        grad = pot%get_gradient('adiabatic', 2)
        d = pot%adiab_dw(:, 2, 2)
        if (size(grad) /= 6) return
        if (any(grad /= [d, 2.0_dp * d, 3.0_dp * d])) return

        ! Coupling vectors are converted the same way.
        ret = 4
        nadv = pot%get_nadv('adiabatic', 1, 2)
        d = pot%adiab_dw(:, 1, 2) / (pot%adiab_w(2, 2) - pot%adiab_w(1, 1))
        if (size(nadv) /= 6) return
        if (any(nadv /= [d, 2.0_dp * d, 3.0_dp * d])) return

        ! The eigenvectors of the previous geometry were saved.
        ret = 5
        if (.not. allocated(pot%ldiab_trans) .or. .not. allocated(pot%lsf_trans)) return
        if (any(pot%ldiab_trans /= trans) .or. any(pot%lsf_trans /= sf_trans)) return
        ret = 0
    end function derivative_test


    integer function group_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        use evaluator_base_mod, only : state_location
        type(group_evaluator) :: pot
        type(state_location) :: loc
        real(dp) :: ones(13, 13), sf_ones(7, 7)

        ones = 1.0_dp
        ret = 1
        ! 4 singlets and 3 triplets, each spin component a separate state.
        call pot%set_spin_manifolds([4, 3], [1, 3], 'components')
        if (pot%n_state /= 13 .or. pot%n_sf /= 7) return
        if (size(pot%manifold) /= 2) return
        if (any(pot%manifold%n /= [4, 3])) return
        if (any(pot%manifold%n_comp /= [1, 3])) return
        if (any(pot%manifold%mult /= [1, 3])) return
        if (any(pot%manifold%i0 /= [1, 5])) return
        if (any(pot%manifold%sf0 /= [1, 5])) return

        ret = 2
        loc = pot%locate_state(4)
        if (loc%group /= 1 .or. loc%manifold /= 1 .or. loc%mult /= 1) return
        if (loc%component /= 1 .or. loc%index /= 4 .or. loc%sf /= 4) return
        loc = pot%locate_state(5)
        if (loc%group /= 2 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 1 .or. loc%index /= 1 .or. loc%sf /= 5) return
        loc = pot%locate_state(10)
        if (loc%group /= 3 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 2 .or. loc%index /= 3 .or. loc%sf /= 7) return
        loc = pot%locate_state(13)
        if (loc%group /= 4 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 3 .or. loc%index /= 3 .or. loc%sf /= 7) return

        ! The spin_free basis has one state for each spin-free state. Summed over the components,
        ! a matrix of ones counts the component pairs.
        ret = 3
        if (pot%basis_size('spin_free') /= 7 .or. pot%basis_size('diabatic') /= 13) return
        sf_ones = pot%sum_components(ones)
        if (sf_ones(1, 4) /= 1.0_dp .or. sf_ones(2, 5) /= 3.0_dp .or. sf_ones(7, 3) /= 3.0_dp) return
        if (sf_ones(5, 7) /= 9.0_dp .or. sum(sf_ones) /= 13.0_dp**2) return
        ret = 0
    end function group_test


    integer function vc_group_test() result(ret)
        use global_defs, only : dp
        use vc_model_mod, only : vc_model
        use vc_evaluator_mod, only : vc_evaluator
        use evaluator_base_mod, only : state_location
        type(vc_model), target :: model
        type(vc_evaluator) :: pot
        type(vc_evaluator) :: pot_bad
        type(state_location) :: loc
        real(dp) :: soc(8, 8)

        ! 2 singlets and 2 triplets, without doublets.
        model%nstate = [2, 0, 2]
        model%nmode = 1
        soc = 0.0_dp
        ! Singlet-triplet couplings only.
        soc(1, 3) = 1.0e-4_dp
        soc(2, 6) = 2.0e-4_dp
        soc(1, 8) = 3.0e-4_dp
        soc = soc + transpose(soc)
        call model%soc%add_zero(model%nmode, soc)

        ret = 1
        call pot%initialize(model)
        if (pot%n_state /= 8 .or. pot%n_sf /= 4) return
        if (any(pot%manifold%n /= [2, 2])) return
        if (any(pot%manifold%n_comp /= [1, 3])) return
        if (any(pot%manifold%i0 /= [1, 3])) return
        if (any(pot%manifold%mult /= [1, 3])) return
        loc = pot%locate_state(6)
        if (loc%manifold /= 2 .or. loc%component /= 2 .or. loc%index /= 2) return

        ret = 2
        if (pot%max_soc_in_groups(model%soc%zero) /= 0.0_dp) return

        ! Coupling between the two states of the second triplet component.
        ret = 3
        soc(5, 6) = 1.0e-4_dp
        soc(6, 5) = 1.0e-4_dp
        call model%soc%add_zero(model%nmode, soc)
        call pot_bad%initialize(model)
        if (pot_bad%max_soc_in_groups(model%soc%zero) <= 0.0_dp) return
        ret = 0
    end function vc_group_test


    !> The spin_free basis of one singlet and two triplets, with couplings h(k, c) between the
    !! singlet and component c of triplet k, and the coupling mix between the two triplets.
    integer function spin_free_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        type(group_evaluator) :: pot
        type(group_evaluator) :: pot_eff
        real(dp) :: h(2, 3), r(3, 3), s2(3, 3), s2_ref(3, 3), u(3, 3)
        real(dp), parameter :: tiny = 1.0e-14_dp
        integer :: k

        h = reshape([1.0e-4_dp, 2.0e-4_dp, -3.0e-4_dp, 0.5e-4_dp, 0.0_dp, 4.0e-4_dp], [2, 3])

        ! Without mixing, the squared coupling of the singlet and triplet k is sum_c h(k, c)**2.
        ret = 1
        call set_st_model(pot, 'components', 0.0_dp, h)
        if (pot%basis_size('spin_free') /= 3 .or. pot%basis_size('adiabatic') /= 7) return
        if (size(pot%get_energy('spin_free')) /= 3) return
        if (any(abs(pot%get_energy('spin_free') - [0.0_dp, 0.1_dp, 0.2_dp]) > tiny)) return
        s2 = pot%get_soc2('spin_free')
        do k = 1, 2
            if (abs(s2(1, k + 1) - sum(h(k, :)**2)) > tiny) return
        end do
        if (any(abs(s2 - transpose(s2)) > tiny)) return

        ! With mixing of the triplets, the couplings are those of the spin-free states, and their
        ! sum over the triplets does not change.
        ret = 2
        call set_st_model(pot, 'components', 0.03_dp, h)
        s2 = pot%get_soc2('spin_free')
        u = pot%get_transformation('diabatic', 'spin_free')
        if (any(shape(u) /= [3, 3]) .or. any(u /= pot%sf_trans)) return
        do k = 1, 2
            if (abs(s2(1, k + 1) - sum(matmul(u(2:3, k + 1), h)**2)) > tiny) return
        end do
        if (abs(s2(1, 2) + s2(1, 3) - sum(h**2)) > tiny) return
        if (abs(pot%get_energy('spin_free', 3) - pot%get_energy('spin_free', 2)                  &
        &       - 2 * sqrt(0.05_dp**2 + 0.03_dp**2)) > tiny) return

        ! A rotation of the triplet components does not change the couplings.
        ret = 3
        s2_ref = s2
        r = reshape([0.6_dp, 0.8_dp, 0.0_dp, -0.8_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp],    &
        &           [3, 3])
        r = matmul(r, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.28_dp, 0.96_dp, 0.0_dp,          &
        &                      -0.96_dp, 0.28_dp], [3, 3]))
        call set_st_model(pot, 'components', 0.03_dp, matmul(h, r))
        if (any(abs(pot%get_soc2('spin_free') - s2_ref) > tiny)) return

        ! Effective states with the root-sum-square couplings give the same couplings without
        ! mixing, and with mixing if the component couplings of the two triplets are parallel.
        ret = 4
        call set_st_model(pot, 'components', 0.0_dp, h)
        call set_st_model(pot_eff, 'effective_rss', 0.0_dp, h)
        if (pot_eff%basis_size('spin_free') /= 3 .or. pot_eff%n_state /= 3) return
        if (any(abs(pot_eff%get_soc2('spin_free') - pot%get_soc2('spin_free')) > tiny)) return
        h(2, :) = 0.5_dp * h(1, :)
        call set_st_model(pot, 'components', 0.03_dp, h)
        call set_st_model(pot_eff, 'effective_rss', 0.03_dp, h)
        if (any(abs(pot_eff%get_soc2('spin_free') - pot%get_soc2('spin_free')) > tiny)) return
        ret = 0
    end function spin_free_test


    !> Set and evaluate the potential of spin_free_test. With effective states, the coupling of
    !! the singlet and triplet k is the root-sum-square of h(k, :).
    subroutine set_st_model(pot, spin_states, mix, h)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        type(group_evaluator), intent(inout) :: pot
        character(len=*), intent(in) :: spin_states
        real(dp), intent(in) :: mix
        real(dp), intent(in) :: h(:, :)
        integer :: k, c, i

        if (.not. allocated(pot%manifold)) then
            call pot%set_spin_manifolds([1, 2], [1, 3], spin_states)
            pot%ndof = 1
            call pot%allocate_hamiltonian()
        end if
        pot%sf_w = 0.0_dp
        pot%sf_w(2, 2) = 0.1_dp
        pot%sf_w(3, 3) = 0.2_dp
        pot%sf_w(2, 3) = mix
        pot%sf_w(3, 2) = mix
        pot%soc_w = 0.0_dp
        do k = 1, 2
            if (spin_states == 'components') then
                do c = 1, 3
                    i = 1 + (c - 1) * 2 + k
                    pot%soc_w(1, i) = h(k, c)
                end do
            else
                pot%soc_w(1, 1 + k) = norm2(h(k, :))
            end if
        end do
        pot%soc_w(:, 1) = pot%soc_w(1, :)
        call pot%eval()
    end subroutine set_st_model


    !> Construction of the effective model, the projection model, the spin-free basis, couplings
    !! of rank 1, the signs of rss and a given reference geometry.
    integer function effective_test() result(ret)
        use global_defs, only : dp
        use matrix_mod, only : unit_mat
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff, eff_p, eff_s, eff_m, eff_sign, eff_ref
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: geom(1, 2), u0(3), c(2, 2), p
        real(dp), allocatable :: e_sf(:), e_src(:)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: n, i, k
        logical :: negative

        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(src, h0, h1, x0, x1)
        geom(1, :) = [0.3_dp, -0.2_dp]

        ! The model has one state for each spin-free state, and P is an isometry.
        ret = 1
        call eff%init(src, 'projection')
        if (eff%n_state /= 4 .or. eff%n_sf /= 4 .or. eff%spin_states /= 'effective') return
        if (any(eff%manifold%n_comp /= [1, 1]) .or. any(eff%manifold%mult /= [1, 3])) return
        if (eff%ndof /= 2 .or. eff%classical_coordinates /= 'model') return
        if (abs(norm2(eff%u) - 1.0_dp) > tol) return
        if (any(shape(eff%proj) /= [8, 4])) return
        if (any(abs(matmul(transpose(eff%proj), eff%proj) - unit_mat(4)) > tol)) return
        if (eff%ref_fraction >= 1.0_dp - 1.0e-6_dp .or. eff%ref_fraction <= 0.0_dp) return

        ! The diabatic Hamiltonian is P^T H P, the spin-free Hamiltonian is that of the source,
        ! and the source builds no bases.
        ret = 2
        call eff%update_geometry(geom, 'model')
        call eff%eval()
        call src%update_geometry(geom, 'model')
        call src%eval()
        if (any(abs(eff%diab_w - matmul(transpose(eff%proj), matmul(src%diab_w, eff%proj)))     &
        &       > tol)) return
        do n = 1, 2
            if (any(abs(eff%diab_dw(n, :, :) - matmul(transpose(eff%proj),                      &
            &       matmul(src%diab_dw(n, :, :), eff%proj))) > tol)) return
        end do
        if (any(eff%sf_w /= src%sf_w) .or. any(eff%sf_dw /= src%sf_dw)) return
        if (any(abs(eff%get_energy('spin_free') - src%get_energy('spin_free')) > tol)) return
        if (allocated(eff%source%adiab_w) .or. allocated(eff%source%sf_trans)) return
        if (abs(eff%min_fraction - min(eff%fraction, 1.0_dp)) > tol) return

        ! Couplings of rank 1: u is their direction, and the adiabatic energies of the source are
        ! those of the model and the dark states at the triplet energies, twice. rss equals
        ! projection, and rms is projection / sqrt(3).
        ret = 3
        u0 = [0.36_dp, 0.48_dp, 0.8_dp]
        c = reshape([1.0e-3_dp, 0.7e-3_dp, -0.5e-3_dp, 0.3e-3_dp], [2, 2])
        do k = 1, 2
            do i = 1, 2
                h0(:, i, k) = c(i, k) * u0
            end do
        end do
        h1 = 0.0_dp
        x0 = 0.0_dp
        x1 = 0.0_dp
        call set_linear_model(src, h0, h1, x0, x1)
        call eff_p%init(src, 'projection')
        call eff_s%init(src, 'rss')
        call eff_m%init(src, 'rms')
        if (any(abs(eff_p%u - u0) > tol) .or. abs(eff_p%ref_fraction - 1.0_dp) > tol) return
        call src%update_geometry(geom, 'model')
        call src%eval()
        e_sf = src%get_energy('spin_free')
        e_src = src%get_energy('adiabatic')
        call eff_p%update_geometry(geom, 'model')
        call eff_p%eval()
        call eff_s%update_geometry(geom, 'model')
        call eff_s%eval()
        call eff_m%update_geometry(geom, 'model')
        call eff_m%eval()
        if (any(abs(sorted([eff_p%get_energy('adiabatic'), e_sf(3:4), e_sf(3:4)]) - e_src)       &
        &       > tol)) return
        if (any(abs(sorted([eff_s%get_energy('adiabatic'), e_sf(3:4), e_sf(3:4)]) - e_src)       &
        &       > tol)) return
        if (any(abs(eff_s%soc_w - eff_p%soc_w) > tol)) return
        if (any(abs(eff_m%soc_w - eff_p%soc_w / sqrt(3.0_dp)) > tol)) return

        ! The rss couplings are |h_ik| with the sign of u.h_ik. Couplings between triplets are
        ! dropped.
        ret = 4
        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(src, h0, h1, x0, x1)
        call eff_sign%init(src, 'rss')
        geom = 0.0_dp
        call eff_sign%update_geometry(geom, 'model')
        call eff_sign%eval()
        negative = .false.
        do k = 1, 2
            do i = 1, 2
                p = dot_product(eff_sign%u, h0(:, i, k))
                negative = negative .or. p < 0.0_dp
                if (abs(eff_sign%soc_w(i, 2 + k) - sign(norm2(h0(:, i, k)), p)) > tol) return
            end do
        end do
        if (.not. negative) return
        if (eff_sign%soc_w(3, 4) /= 0.0_dp) return

        ! Without couplings at the default reference (x = 0), u is computed at the given one.
        ret = 5
        h0 = 0.0_dp
        h1 = 0.0_dp
        do k = 1, 2
            do i = 1, 2
                h1(:, i, k, 1) = c(i, k) * [0.6_dp, 0.0_dp, -0.8_dp]
            end do
        end do
        call set_linear_model(src, h0, h1, x0, x1)
        call eff_ref%init(src, 'projection', reshape([1.0_dp, 0.0_dp], [1, 2]), 'model')
        if (any(abs(eff_ref%u - [-0.6_dp, 0.0_dp, 0.8_dp]) > tol)) return
        if (abs(eff_ref%ref_fraction - 1.0_dp) > tol) return
        ret = 0
    end function effective_test


    !> A rotation of the triplet components rotates u with them and does not change the model.
    integer function effective_rotation_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src, src_rot
        type(effective_spin_evaluator) :: eff(3), eff_rot(3)
        character(len=10), parameter :: models(3) = [character(len=10) :: 'projection', 'rss', 'rms']
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: h0_rot(3, 2, 2), h1_rot(3, 2, 2, 2), r(3, 3), geom(1, 2)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: m, i, k, n

        call general_couplings(h0, h1, x0, x1)
        ! Couplings between triplets would couple the states within a group after the rotation.
        x0 = 0.0_dp
        x1 = 0.0_dp
        r = reshape([0.6_dp, 0.8_dp, 0.0_dp, -0.8_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp],    &
        &           [3, 3])
        r = matmul(r, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.28_dp, 0.96_dp, 0.0_dp,          &
        &                      -0.96_dp, 0.28_dp], [3, 3]))
        do k = 1, 2
            do i = 1, 2
                h0_rot(:, i, k) = matmul(r, h0(:, i, k))
                do n = 1, 2
                    h1_rot(:, i, k, n) = matmul(r, h1(:, i, k, n))
                end do
            end do
        end do
        call set_linear_model(src, h0, h1, x0, x1)
        call set_linear_model(src_rot, h0_rot, h1_rot, x0, x1)
        geom(1, :) = [0.3_dp, -0.2_dp]

        do m = 1, 3
            ret = m
            call eff(m)%init(src, trim(models(m)))
            call eff_rot(m)%init(src_rot, trim(models(m)))
            if (abs(abs(dot_product(eff_rot(m)%u, matmul(r, eff(m)%u))) - 1.0_dp) > tol) return
            call eff(m)%update_geometry(geom, 'model')
            call eff(m)%eval()
            call eff_rot(m)%update_geometry(geom, 'model')
            call eff_rot(m)%eval()
            if (any(abs(eff(m)%get_energy('adiabatic') - eff_rot(m)%get_energy('adiabatic'))     &
            &       > tol)) return
            if (any(abs(abs(eff(m)%soc_w) - abs(eff_rot(m)%soc_w)) > tol)) return
        end do
        ret = 0
    end function effective_rotation_test


    !> Gradients of the adiabatic states of each model match finite differences.
    integer function effective_gradient_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff(6), eff_d
        character(len=10), parameter :: models(6) = [character(len=10) :: 'projection', 'rss',    &
        &                                            'rms', 'sq3', 'rss', 'sq3']
        logical, parameter :: signed(6) = [.true., .true., .true., .true., .false., .false.]
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: geom(1, 2), geom_d(1, 2), e_p(4), e_m(4)
        real(dp), allocatable :: grad(:)
        real(dp), parameter :: delta = 1.0e-5_dp, tol = 1.0e-8_dp
        integer :: m, n, s

        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(src, h0, h1, x0, x1)
        geom(1, :) = [0.3_dp, -0.2_dp]

        do m = 1, 6
            ret = m
            call eff(m)%init(src, trim(models(m)), signed=signed(m))
            call eff(m)%update_geometry(geom, 'model')
            call eff(m)%eval()
            do n = 1, 2
                eff_d = eff(m)
                geom_d = geom
                geom_d(1, n) = geom(1, n) + delta
                call eff_d%update_geometry(geom_d, 'model')
                call eff_d%eval()
                e_p = eff_d%get_energy('adiabatic')
                geom_d(1, n) = geom(1, n) - delta
                call eff_d%update_geometry(geom_d, 'model')
                call eff_d%eval()
                e_m = eff_d%get_energy('adiabatic')
                do s = 1, 4
                    grad = eff(m)%get_gradient('adiabatic', s)
                    if (abs(grad(n) - (e_p(s) - e_m(s)) / (2.0_dp * delta)) > tol) return
                end do
            end do
        end do
        ret = 0
    end function effective_gradient_test


    !----------------------------------------------------------------------------------------------
    !> Signs of the rss, rms and sq3 couplings by continuity, and the unsigned variants.
    !! Along x = 0, 0.1, ..., 2:
    !! - h_11 = (1 - x) a + x c, with c perpendicular to a, rotates past the plane u.h = 0 at
    !!   x = 1 without vanishing. Its coupling keeps the sign of the first point.
    !! - h_21 = (0.5 - x) d, with d parallel to a, is exactly zero at x = 0.5. Its coupling is
    !!   (0.5 - x) |d| when signed, and |0.5 - x| |d| without signs.
    !! A second evaluation at the same geometry gives the same couplings.
    !----------------------------------------------------------------------------------------------
    integer function effective_sign_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff(5)
        character(len=3), parameter :: models(5) = ['rss', 'rms', 'sq3', 'rss', 'sq3']
        logical, parameter :: signed(5) = [.true., .true., .true., .false., .false.]
        real(dp), parameter :: scale(5) = [1.0_dp, 1.0_dp / sqrt(3.0_dp), 1.0_dp / 3.0_dp,        &
        &                                  1.0_dp, 1.0_dp / 3.0_dp]
        real(dp), parameter :: a(3) = [0.6e-3_dp, 0.0_dp, 0.8e-3_dp], c(3) = [0.0_dp, 1.0e-3_dp,  &
        &                      0.0_dp], d(3) = [1.2e-3_dp, 0.0_dp, 1.6e-3_dp]
        real(dp), parameter :: tol = 1.0e-15_dp
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: geom(1, 2), x, c11, c21, dc21
        real(dp), allocatable :: soc_w(:, :)
        integer :: m, i

        h0 = 0.0_dp
        h1 = 0.0_dp
        x0 = 0.0_dp
        x1 = 0.0_dp
        h0(:, 1, 1) = a
        h1(:, 1, 1, 1) = c - a
        h0(:, 2, 1) = 0.5_dp * d
        h1(:, 2, 1, 1) = -d
        call set_linear_model(src, h0, h1, x0, x1)

        do m = 1, 5
            ret = m
            call eff(m)%init(src, models(m), signed=signed(m))
            if (any(abs(eff(m)%u - a / norm2(a)) > tol * 1.0e3_dp)) return
            do i = 0, 20
                x = real(i, dp) / 10.0_dp
                geom(1, :) = [x, 0.0_dp]
                call eff(m)%update_geometry(geom, 'model')
                call eff(m)%eval()
                c11 = scale(m) * norm2((1.0_dp - x) * a + x * c)
                c21 = scale(m) * (0.5_dp - x) * norm2(d)
                dc21 = -scale(m) * norm2(d)
                if (.not. signed(m)) then
                    c21 = abs(c21)
                    dc21 = sign(dc21, x - 0.5_dp)
                end if
                if (abs(eff(m)%soc_w(1, 3) - c11) > tol) return
                if (abs(eff(m)%soc_w(2, 3) - c21) > tol) return
                if (i /= 5 .and. abs(eff(m)%soc_dw(1, 2, 3) - dc21) > tol) return
                if (any(eff(m)%soc_w(:, 4) /= 0.0_dp)) return
                ! A second evaluation, with and without saving the current point.
                if (i == 5 .or. i == 6) then
                    soc_w = eff(m)%soc_w
                    call eff(m)%eval()
                    if (any(eff(m)%soc_w /= soc_w)) return
                    call eff(m)%update_geometry(geom, 'model')
                    call eff(m)%eval()
                    if (any(eff(m)%soc_w /= soc_w)) return
                end if
            end do
        end do
        ret = 0
    end function effective_sign_test


    !> A copy of the wrapper is independent of the original, and derivatives are converted to the
    !! target coordinates by the source.
    integer function effective_copy_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator, pot_holder
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff, eff_scaled
        type(pot_holder) :: a, b
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: geom(1, 2)
        real(dp), allocatable :: grad(:), nadv(:), d(:)

        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(src, h0, h1, x0, x1)
        call eff%init(src, 'projection')
        geom(1, :) = [0.3_dp, -0.2_dp]
        call eff%update_geometry(geom, 'model')
        call eff%eval()

        ! Copy as in the trajectory data, and move the copy to another geometry.
        ret = 1
        allocate(a%pot, source=eff)
        b = a
        select type(p => b%pot)
        type is (effective_spin_evaluator)
            geom(1, :) = [-0.4_dp, 0.1_dp]
            call p%update_geometry(geom, 'model')
            call p%eval()
            if (all(abs(p%adiab_w - eff%adiab_w) < 1.0e-12_dp)) return
        class default
            return
        end select
        select type(p => a%pot)
        type is (effective_spin_evaluator)
            if (any(p%adiab_w /= eff%adiab_w)) return
            if (any(p%source%sf_w /= eff%source%sf_w)) return
            if (any(p%source%soc_w /= eff%source%soc_w)) return
            select type(s => p%source)
            type is (linear_evaluator)
                if (any(s%x /= [0.3_dp, -0.2_dp])) return
            class default
                return
            end select
        class default
            return
        end select

        ! Derivatives are converted by the source.
        ret = 2
        src%scaled = .true.
        call eff_scaled%init(src, 'projection')
        call eff_scaled%update_geometry(geom, 'model')
        call eff_scaled%eval()
        grad = eff_scaled%get_gradient('adiabatic', 1)
        d = eff_scaled%adiab_dw(:, 1, 1)
        if (size(grad) /= 6) return
        if (any(grad /= [d, 2.0_dp * d, 3.0_dp * d])) return
        nadv = eff_scaled%get_nadv('adiabatic', 1, 2)
        d = eff_scaled%adiab_dw(:, 1, 2) / (eff_scaled%adiab_w(2, 2) - eff_scaled%adiab_w(1, 1))
        if (size(nadv) /= 6) return
        if (any(nadv /= [d, 2.0_dp * d, 3.0_dp * d])) return
        ret = 0
    end function effective_copy_test


    !> The Granucci models: E_S and E_T cross at r_c = 10, the couplings change sign between
    !! r_sigma -+ dr_sigma / 2 with the strength of the article, the angle between the couplings
    !! of the two singlets is phi, and the derivatives match finite differences. The retained
    !! fraction of the effective model of granucci-sst is (1 + cos phi) / 2.
    integer function model_1d_test() result(ret)
        use global_defs, only : dp, pi
        use string_mod, only : string
        use model_1d_mod, only : model_1d_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(model_1d_evaluator) :: st, sst
        type(effective_spin_evaluator) :: eff
        real(dp), parameter :: tol = 1.0e-15_dp
        real(dp), parameter :: gamma = sqrt(2.0_dp) * 1.0e-3_dp
        real(dp) :: h(3), h1(3), h2(3)

        ret = 1
        call st%init([string('granucci-st'), string('r_sigma'), string('9.5')])
        if (st%n_state /= 4 .or. st%n_sf /= 2 .or. any(st%manifold%mult /= [1, 3])) return
        if (st%spin_states /= 'components') return

        ! The energies cross at r_c, and a1 is that of the article.
        call eval_model(st, 10.0_dp)
        if (abs(st%sf_w(1, 1) - st%sf_w(2, 2)) > tol) return
        if (abs(st%sf_w(1, 1) - 0.04_dp - 0.03452_dp * exp(-3.5_dp)) > 1.0e-7_dp) return

        ! The couplings are h before r_sigma - dr_sigma / 2, zero at r_sigma and -h after
        ! r_sigma + dr_sigma / 2.
        ret = 2
        call eval_model(st, 8.5_dp)
        h = st%soc_w(1, 2:4)
        if (abs(norm2(h) - gamma) > tol) return
        if (any(st%soc_w(2:4, 1) /= h) .or. any(st%soc_w(2:4, 2:4) /= 0.0_dp)) return
        call eval_model(st, 9.5_dp)
        if (any(abs(st%soc_w(1, 2:4)) > tol)) return
        call eval_model(st, 10.5_dp)
        if (any(st%soc_w(1, 2:4) /= -h)) return

        ! Derivatives, inside and outside the sign change.
        ret = 3
        if (.not. derivatives_ok(st, [8.0_dp, 9.2_dp, 9.9_dp, 12.0_dp])) return

        ! Two singlets with couplings of equal strength at the angle phi.
        ret = 4
        call sst%init([string('granucci-sst'), string('PHI'), string('60')])
        if (sst%n_state /= 5 .or. sst%n_sf /= 3 .or. any(sst%manifold%n /= [2, 1])) return
        call eval_model(sst, 11.0_dp)
        h1 = sst%soc_w(1, 3:5)
        h2 = sst%soc_w(2, 3:5)
        if (abs(norm2(h1) - gamma) > tol .or. abs(norm2(h2) - gamma) > tol) return
        if (abs(dot_product(h1, h2) / gamma**2 - 0.5_dp) > 1.0e-14_dp) return
        if (abs(sst%sf_w(2, 2) - sst%sf_w(1, 1) - 0.035_dp) > tol) return
        if (sst%sf_w(1, 2) /= 0.0_dp) return
        if (.not. derivatives_ok(sst, [9.0_dp, 15.0_dp])) return

        ! The retained fraction of the effective model.
        ret = 5
        call eff%init(sst, 'projection')
        if (abs(eff%ref_fraction - (1.0_dp + cos(pi / 3.0_dp)) / 2.0_dp) > 1.0e-14_dp) return
        ret = 0
    end function model_1d_test


    !> Evaluate a model_1d_evaluator at position x.
    subroutine eval_model(model, x)
        use global_defs, only : dp
        use model_1d_mod, only : model_1d_evaluator
        type(model_1d_evaluator), intent(inout) :: model
        real(dp), intent(in) :: x

        call model%update_geometry(reshape([x], [1, 1]), 'model')
        call model%eval()
    end subroutine eval_model


    !> Central finite differences of sf_w and soc_w of a model_1d_evaluator match sf_dw and
    !! soc_dw at each point.
    logical function derivatives_ok(model, points) result(ok)
        use global_defs, only : dp
        use model_1d_mod, only : model_1d_evaluator
        type(model_1d_evaluator), intent(inout) :: model
        real(dp), intent(in) :: points(:)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp), allocatable :: sf_p(:, :), soc_p(:, :)
        integer :: i

        ok = .false.
        do i = 1, size(points)
            call eval_model(model, points(i) + step)
            sf_p = model%sf_w
            soc_p = model%soc_w
            call eval_model(model, points(i) - step)
            sf_p = (sf_p - model%sf_w) / (2.0_dp * step)
            soc_p = (soc_p - model%soc_w) / (2.0_dp * step)
            call eval_model(model, points(i))
            if (any(abs(sf_p - model%sf_dw(1, :, :)) > 1.0e-10_dp)) return
            if (any(abs(soc_p - model%soc_dw(1, :, :)) > 1.0e-10_dp)) return
        end do
        ok = .true.
    end function derivatives_ok


    !> The state layout of a potential with components, of a potential with effective states and
    !! of the effective model built from the first, which adds the model and the layout of its
    !! source.
    integer function state_layout_test() result(ret)
        use string_mod, only : string
        use json_module, only : json_core, json_value
        use pot_test_mod, only : group_evaluator
        use model_1d_mod, only : model_1d_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(model_1d_evaluator) :: st
        type(group_evaluator) :: grp
        type(effective_spin_evaluator) :: eff
        type(json_core) :: json
        type(json_value), pointer :: root

        call json%initialize()

        ! Components: 1 singlet and the 3 components of 1 triplet.
        ret = 1
        call st%init([string('granucci-st')])
        call json%create_object(root, '')
        call st%add_state_layout(json, root)
        if (jstr(json, root, 'spin_states') /= 'components') return
        if (jint(json, root, 'n_state.diabatic') /= 4) return
        if (jint(json, root, 'n_state.adiabatic') /= 4) return
        if (jint(json, root, 'n_state.spin_free') /= 2) return
        if (jint(json, root, 'manifolds(1).multiplicity') /= 1) return
        if (jint(json, root, 'manifolds(1).n') /= 1) return
        if (jint(json, root, 'manifolds(2).multiplicity') /= 3) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 3) return
        if (jint(json, root, 'manifolds(2).first_spin_free') /= 2) return
        if (jint(json, root, 'manifolds(2).first_state') /= 2) return
        call json%destroy(root)

        ! Effective states given as input.
        ret = 2
        call grp%set_spin_manifolds([1, 2], [1, 3], 'effective_rss')
        call json%create_object(root, '')
        call grp%add_state_layout(json, root)
        if (jstr(json, root, 'spin_states') /= 'effective') return
        if (jint(json, root, 'n_state.diabatic') /= 3) return
        if (jint(json, root, 'n_state.spin_free') /= 3) return
        if (jint(json, root, 'manifolds(2).n') /= 2) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 1) return
        call json%destroy(root)

        ! Effective model: its own layout, the model, and the layout of the source.
        ret = 3
        call eff%init(st, 'sq3', signed=.false.)
        call json%create_object(root, '')
        call eff%add_state_layout(json, root)
        if (jstr(json, root, 'spin_states') /= 'effective') return
        if (jint(json, root, 'n_state.diabatic') /= 2) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 1) return
        if (jstr(json, root, 'effective_spin_states.model') /= 'sq3') return
        if (jlog(json, root, 'effective_spin_states.signed')) return
        if (jstr(json, root, 'source.spin_states') /= 'components') return
        if (jint(json, root, 'source.n_state.diabatic') /= 4) return
        if (jint(json, root, 'source.manifolds(2).n_comp') /= 3) return
        call json%destroy(root)
        if (json%failed()) return
        ret = 0
    end function state_layout_test


    !> Integer at path in a JSON object, or -1 if not found.
    integer function jint(json, root, path) result(val)
        use json_module, only : json_core, json_value
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: root
        character(len=*), intent(in) :: path
        logical :: found

        call json%get(root, path, val, found)
        if (.not. found) val = -1
    end function jint


    !> Logical at path in a JSON object, or .true. if not found.
    logical function jlog(json, root, path) result(val)
        use json_module, only : json_core, json_value
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: root
        character(len=*), intent(in) :: path
        logical :: found

        call json%get(root, path, val, found)
        if (.not. found) val = .true.
    end function jlog


    !> String at path in a JSON object, or an empty string if not found.
    function jstr(json, root, path) result(val)
        use json_module, only : json_core, json_value
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: root
        character(len=*), intent(in) :: path
        character(len=:), allocatable :: val
        logical :: found

        call json%get(root, path, val, found)
        if (.not. found) val = ''
    end function jstr


    !> Cases which stop with an error: in effective_init, a doublet, a source of effective states,
    !! no singlet-triplet coupling at the reference geometry and the projection model without
    !! signs; in model_init, an unknown parameter.
    subroutine effective_error(case)
        use global_defs, only : dp
        use string_mod, only : string
        use pot_test_mod, only : linear_evaluator
        use model_1d_mod, only : model_1d_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        character(len=*), intent(in) :: case
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff
        type(model_1d_evaluator) :: model
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)

        select case(case)
        case('model_parameter')
            call model%init([string('granucci-st'), string('r_sigma'), string('9.5'),            &
            &                string('sigma'), string('1.0')])
            return
        case('doublet')
            call src%set_spin_manifolds([1, 1], [1, 2], 'components')
            src%classical_coordinates = 'model'
        case('effective')
            call src%set_spin_manifolds([1, 1], [1, 3], 'effective_rss')
            src%classical_coordinates = 'model'
        case('no_coupling')
            h0 = 0.0_dp
            h1 = 0.0_dp
            h1(1, 1, 1, 1) = 1.0e-4_dp
            x0 = 0.0_dp
            x1 = 0.0_dp
            call set_linear_model(src, h0, h1, x0, x1)
        case('unsigned_projection')
            call general_couplings(h0, h1, x0, x1)
            call set_linear_model(src, h0, h1, x0, x1)
            call eff%init(src, 'projection', signed=.false.)
            return
        case default
            write(*, *) 'Unknown case: ', case
            return
        end select
        call eff%init(src, 'projection')
    end subroutine effective_error


    !> Couplings h0(c, i, k) of singlet i and component c of triplet k, their derivatives
    !! h1(c, i, k, n), and the couplings x0(c, c') of component c of T1 and component c' of T2
    !! and their derivatives x1(c, c', n). Not of rank 1, and zero within a group.
    subroutine general_couplings(h0, h1, x0, x1)
        use global_defs, only : dp
        real(dp), intent(out) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        integer :: c, cc, i, k, n

        do k = 1, 2
            do i = 1, 2
                do c = 1, 3
                    h0(c, i, k) = 1.0e-3_dp * cos(real(c + 2 * i + 4 * k, dp))
                    do n = 1, 2
                        h1(c, i, k, n) = 2.0e-4_dp * sin(real(c + 3 * i + 5 * k + 7 * n, dp))
                    end do
                end do
            end do
        end do
        x0 = 0.0_dp
        x1 = 0.0_dp
        do cc = 1, 3
            do c = 1, 3
                if (c == cc) cycle
                x0(c, cc) = 3.0e-4_dp * sin(real(c + 2 * cc, dp))
                do n = 1, 2
                    x1(c, cc, n) = 1.0e-4_dp * cos(real(c + cc + n, dp))
                end do
            end do
        end do
    end subroutine general_couplings


    !> Linear model of two singlets and two triplets in two coordinates, with the couplings of
    !! general_couplings. States of the source: S1, S2, T1(1), T2(1), T1(2), T2(2), T1(3), T2(3).
    subroutine set_linear_model(pot, h0, h1, x0, x1)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        type(linear_evaluator), intent(out) :: pot
        real(dp), intent(in) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp), parameter :: slope(4) = [0.02_dp, -0.03_dp, 0.01_dp, 0.015_dp]
        integer :: i, n

        call pot%set_spin_manifolds([2, 2], [1, 3], 'components')
        pot%ndof = 2
        pot%classical_coordinates = 'model'
        pot%x = [0.0_dp, 0.0_dp]
        pot%sf0 = reshape([0.0_dp, 0.01_dp, 0.0_dp, 0.0_dp,                                       &
        &                  0.01_dp, 0.1_dp, 0.0_dp, 0.0_dp,                                       &
        &                  0.0_dp, 0.0_dp, 0.06_dp, 0.005_dp,                                     &
        &                  0.0_dp, 0.0_dp, 0.005_dp, 0.09_dp], [4, 4])
        allocate(pot%sf1(2, 4, 4), source=0.0_dp)
        do i = 1, 4
            pot%sf1(1, i, i) = slope(i)
        end do
        pot%sf1(2, 1, 2) = 0.004_dp
        pot%sf1(2, 2, 1) = 0.004_dp
        pot%sf1(2, 3, 4) = 0.002_dp
        pot%sf1(2, 4, 3) = 0.002_dp
        pot%sf1(2, 3, 3) = -0.01_dp
        pot%soc0 = st_matrix(h0, x0)
        allocate(pot%soc1(2, 8, 8))
        do n = 1, 2
            pot%soc1(n, :, :) = st_matrix(h1(:, :, :, n), x1(:, :, n))
        end do
    end subroutine set_linear_model


    !> Spin-orbit couplings of the linear model, from the singlet-triplet couplings h(c, i, k) and
    !! the couplings x(c, c') of component c of T1 and component c' of T2.
    pure function st_matrix(h, x) result(soc)
        use global_defs, only : dp
        real(dp), intent(in) :: h(3, 2, 2), x(3, 3)
        real(dp) :: soc(8, 8)
        integer :: i, k, c, cc

        soc = 0.0_dp
        do k = 1, 2
            do c = 1, 3
                do i = 1, 2
                    soc(i, 2 + (c - 1) * 2 + k) = h(c, i, k)
                end do
            end do
        end do
        do cc = 1, 3
            do c = 1, 3
                soc(2 + (c - 1) * 2 + 1, 2 + (cc - 1) * 2 + 2) = x(c, cc)
            end do
        end do
        soc = soc + transpose(soc)
    end function st_matrix


    pure function sorted(v) result(s)
        use global_defs, only : dp
        real(dp), intent(in) :: v(:)
        real(dp) :: s(size(v))
        real(dp) :: t
        integer :: i, j

        s = v
        do i = 2, size(s)
            t = s(i)
            j = i - 1
            do while (j >= 1)
                if (s(j) <= t) exit
                s(j + 1) = s(j)
                j = j - 1
            end do
            s(j + 1) = t
        end do
    end function sorted


end program pot_test
