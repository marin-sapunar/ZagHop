module pot_test_mod
    use global_defs, only : dp
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

contains

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


end program pot_test
