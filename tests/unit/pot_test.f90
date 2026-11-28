module pot_test_mod
    use global_defs, only : dp
    use evaluator_base_mod, only : potential_evaluator
    use evaluator_diabatic_mod, only : diabatic_evaluator
    implicit none

    !> Evaluator without a potential, used to test the spin manifold setup.
    type, extends(diabatic_evaluator) :: group_evaluator
    contains


        procedure :: update_geometry => group_update_geometry
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
        procedure :: update_geometry => linear_update_geometry
        procedure :: eval_diab => linear_eval_diab
        procedure :: default_reference => linear_default_reference
        procedure :: get_derivative_target_coord => linear_derivative
    end type linear_evaluator

    !> Potential stored as in the trajectory data, to test copies.
    type :: pot_holder
        class(potential_evaluator), allocatable :: pot
    end type pot_holder

contains

    subroutine linear_update_geometry(self, geometry, coordinates)
        class(linear_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        self%x = geometry(1, :)
    end subroutine linear_update_geometry

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

    subroutine group_update_geometry(self, geometry, coordinates)
        class(group_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
    end subroutine group_update_geometry

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
        if (index(arg, 'external_') == 1) then
            call external_error(trim(arg))
        else if (index(arg, 'transform_') == 1) then
            call transform_error(trim(arg))
        else
            call effective_error(trim(arg))
        end if
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

    ret = mch_test()
    if (ret /= 0) then
        write(*,*) "Error in mch_test."
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

    ret = link_test()
    if (ret /= 0) then
        write(*,*) "Error in link_test."
        stop ret
    end if

    ret = layout_subset_test()
    if (ret /= 0) then
        write(*,*) "Error in layout_subset_test."
        stop ret
    end if

    ret = subset_test()
    if (ret /= 0) then
        write(*,*) "Error in subset_test."
        stop ret
    end if

    ret = subset_gradient_test()
    if (ret /= 0) then
        write(*,*) "Error in subset_gradient_test."
        stop ret
    end if

    ret = effective_subset_test()
    if (ret /= 0) then
        write(*,*) "Error in effective_subset_test."
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

    ret = save_state_test()
    if (ret /= 0) then
        write(*,*) "Error in save_state_test."
        stop ret
    end if

    ret = external_json_test()
    if (ret /= 0) then
        write(*,*) "Error in external_json_test."
        stop ret
    end if

    ret = basis_test()
    if (ret /= 0) then
        write(*,*) "Error in basis_test."
        stop ret
    end if

    ret = populations_test()
    if (ret /= 0) then
        write(*,*) "Error in populations_test."
        stop ret
    end if

    ret = transform_test()
    if (ret /= 0) then
        write(*,*) "Error in transform_test."
        stop ret
    end if

    ! A test which stops in the tested code ends the program with exit status 0. This line is
    ! required by tests/unit/run_cases.cmake.
    write(*, *) 'All tests passed.'

    contains


    !> Populations from the coefficients, |transform(from, to, c)|**2.
    !! - Singlet and two triplets of mch_test, with spin-orbit couplings: from normalised
    !!   adiabatic coefficients c, the populations of the adiabatic, mch and diabatic bases sum to
    !!   1, the adiabatic ones are |c|**2 and the diabatic ones |U c|**2. From the diabatic
    !!   coefficients U c, the populations are the same.
    !! - The mch populations summed over the spin components (sum_components) are
    !!   sum_c |sum_k mch_trans(k, a) d(k, c)|**2, with d the diabatic coefficients.
    !! - The same without spin-orbit couplings: the adiabatic states are the MCH states, so the
    !!   summed populations of both bases are the same for each MCH state.
    !! - A rotation of the triplet components of the diabatic coefficients changes the mch
    !!   populations of the components, but not their sums.
    !! - S1 and T1 of the linear model (subset): the populations of all bases sum to 1.
    !! - Effective spin states: the mch populations are those of the spin-free states.
    integer function populations_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator, linear_evaluator
        type(group_evaluator) :: pot, pot0, pot_eff
        type(linear_evaluator) :: sub
        real(dp) :: h(2, 3), r(3, 3), h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp), allocatable :: u(:, :), p_ad(:), p_di(:), p_mch(:), p_sf(:), p_ref(:)
        real(dp), allocatable :: e_ad(:), e_mch(:), sum_ad(:)
        complex(dp), allocatable :: c(:), d(:), d_rot(:)
        complex(dp) :: amp
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: i, k, cc, a

        h = reshape([1.0e-4_dp, 2.0e-4_dp, -3.0e-4_dp, 0.5e-4_dp, 0.0_dp, 4.0e-4_dp], [2, 3])

        ! Sums, adiabatic and diabatic populations.
        ret = 1
        call set_st_model(pot, 'components', 0.03_dp, h)
        c = test_coefficients(7)
        u = pot%get_transformation('diabatic', 'adiabatic')
        p_ad = abs(pot%transform('adiabatic', 'adiabatic', c))**2
        p_di = abs(pot%transform('adiabatic', 'diabatic', c))**2
        p_mch = abs(pot%transform('adiabatic', 'mch', c))**2
        if (size(p_ad) /= 7 .or. size(p_di) /= 7 .or. size(p_mch) /= 7) return
        if (abs(sum(p_ad) - 1.0_dp) > tol .or. abs(sum(p_di) - 1.0_dp) > tol) return
        if (abs(sum(p_mch) - 1.0_dp) > tol) return
        if (any(abs(p_ad - abs(c)**2) > tol)) return
        d = matmul(u, c)
        if (any(abs(p_di - abs(d)**2) > tol)) return
        ret = 2
        if (any(abs(abs(pot%transform('diabatic', 'adiabatic', d))**2 - p_ad) > tol)) return
        if (any(abs(abs(pot%transform('diabatic', 'diabatic', d))**2 - p_di) > tol)) return
        if (any(abs(abs(pot%transform('diabatic', 'mch', d))**2 - p_mch) > tol)) return

        ! Summed over the components (state 1 + (c - 1) * 2 + k is component c of T_k).
        ret = 3
        p_sf = pot%states%sum_components(p_mch)
        allocate(p_ref(3), source=0.0_dp)
        p_ref(1) = abs(d(1))**2
        do a = 2, 3
            do cc = 1, 3
                amp = sum(pot%mch_trans(2:3, a) * d(1 + (cc - 1) * 2 + [1, 2]))
                p_ref(a) = p_ref(a) + abs(amp)**2
            end do
        end do
        if (size(p_sf) /= 3 .or. any(abs(p_sf - p_ref) > tol)) return
        if (abs(sum(p_sf) - 1.0_dp) > tol) return

        ! Without spin-orbit couplings, the adiabatic states are the MCH states.
        ret = 4
        call set_st_model(pot0, 'components', 0.03_dp, 0.0_dp * h)
        p_ad = abs(pot0%transform('adiabatic', 'adiabatic', c))**2
        p_sf = pot0%states%sum_components(abs(pot0%transform('adiabatic', 'mch', c))**2)
        e_ad = pot0%get_energy('adiabatic')
        e_mch = pot0%get_energy('mch_sf')
        allocate(sum_ad(3), source=0.0_dp)
        do i = 1, 7
            a = minloc(abs(e_mch - e_ad(i)), dim=1)
            if (abs(e_mch(a) - e_ad(i)) > 1.0e-10_dp) return
            sum_ad(a) = sum_ad(a) + p_ad(i)
        end do
        if (any(abs(p_sf - sum_ad) > tol)) return

        ! Rotation of the triplet components.
        ret = 5
        r = reshape([0.6_dp, 0.8_dp, 0.0_dp, -0.8_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp],    &
        &           [3, 3])
        r = matmul(r, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.28_dp, 0.96_dp, 0.0_dp,          &
        &                      -0.96_dp, 0.28_dp], [3, 3]))
        d = test_coefficients(7)
        d_rot = d
        do k = 1, 2
            do cc = 1, 3
                d_rot(1 + (cc - 1) * 2 + k) = sum(r(cc, :) * d(1 + [0, 2, 4] + k))
            end do
        end do
        if (abs(sum(abs(d_rot)**2) - 1.0_dp) > tol) return
        p_mch = abs(pot%transform('diabatic', 'mch', d))**2
        p_ref = abs(pot%transform('diabatic', 'mch', d_rot))**2
        if (all(abs(p_ref - p_mch) < 1.0e-3_dp)) return
        if (any(abs(pot%states%sum_components(p_ref) - pot%states%sum_components(p_mch))         &
        &       > tol)) return

        ! Subset: 4 adiabatic states, 8 diabatic states, 4 MCH states (2 multiplets).
        ret = 6
        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(sub, h0, h1, x0, x1)
        call sub%select_states([1, 1], [1, 3])
        call sub%update_geometry(reshape([0.1_dp, 0.2_dp], [1, 2]), 'model')
        call sub%eval()
        c = test_coefficients(4)
        p_di = abs(sub%transform('adiabatic', 'diabatic', c))**2
        p_mch = abs(sub%transform('adiabatic', 'mch', c))**2
        if (size(p_di) /= 8 .or. size(p_mch) /= 4) return
        if (abs(sum(p_di) - 1.0_dp) > tol .or. abs(sum(p_mch) - 1.0_dp) > tol) return
        if (size(sub%states%sum_components(p_mch)) /= 2) return
        d = matmul(sub%get_transformation('diabatic', 'adiabatic'), c)
        if (any(abs(abs(sub%transform('diabatic', 'adiabatic', d))**2 - abs(c)**2) > tol)) return

        ! Effective spin states: one state for each multiplet in both spaces.
        ret = 7
        call set_st_model(pot_eff, 'effective_rss', 0.03_dp, h)
        c = test_coefficients(3)
        p_mch = abs(pot_eff%transform('adiabatic', 'mch', c))**2
        if (size(p_mch) /= 3 .or. abs(sum(p_mch) - 1.0_dp) > tol) return
        if (any(pot_eff%states%sum_components(p_mch) /= p_mch)) return
        ret = 0
    end function populations_test


    !----------------------------------------------------------------------------------------------
    !> Transformations between the bases (can_transform, get_transformation, transform).
    !! - can_transform for the singlet and two triplets of mch_test (all states selected), the
    !!   linear model with S1 and T1 selected, effective spin states and an external program:
    !!   only within a space, and without the diabatic bases for the external program.
    !!   get_transformation gives a matrix with the sizes of the two bases for each pair where
    !!   can_transform is true.
    !! - Composition: diabatic -> adiabatic is adiab_trans of the evaluator and the product of
    !!   diabatic -> mch and mch -> adiabatic; the columns of each are orthonormal, also for the
    !!   subset (rectangular); diabatic_sf -> mch_sf is the same for every component of diabatic
    !!   -> mch.
    !! - Coefficients: transform gives T^T c, also for the subset. Coefficients in the adiabatic
    !!   states, transformed to the diabatic or mch basis and back, are unchanged.
    !! - Operators: the Hamiltonians of the diabatic and mch bases are diagonal in the adiabatic
    !!   basis, with the adiabatic energies (U^T H U = V^T H_sub V), and the spin-free Hamiltonian
    !!   is diagonal in the mch_sf basis. The Hamiltonian of the mch basis has the energies of the
    !!   spin-free states on the diagonal.
    !! - Without spin-orbit couplings, mch -> adiabatic takes each adiabatic state to one MCH
    !!   state of the same energy.
    !! - Equal bases: both specifics return their argument. Effective spin states: the
    !!   transformations of both spaces are the same.
    !----------------------------------------------------------------------------------------------
    integer function transform_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator, linear_evaluator
        use external_json_mod, only : external_json_evaluator
        use matrix_mod, only : diag, unit_mat
        type(group_evaluator) :: pot, pot0, pot_eff
        type(linear_evaluator) :: sub
        type(external_json_evaluator) :: ext
        !> can_transform(bases(i), bases(j)) for the bases diabatic, mch, adiabatic, diabatic_sf,
        !! mch_sf: for potentials with diabatic states, and for the external evaluator.
        logical, parameter :: y = .true., n = .false.
        logical, parameter :: diab_table(5, 5) = reshape([y, y, y, n, n,                         &
        &                                                 y, y, y, n, n,                         &
        &                                                 y, y, y, n, n,                         &
        &                                                 n, n, n, y, y,                         &
        &                                                 n, n, n, y, y], [5, 5])
        logical, parameter :: ext_table(5, 5) = reshape([n, n, n, n, n,                          &
        &                                                n, y, y, n, n,                          &
        &                                                n, y, y, n, n,                          &
        &                                                n, n, n, n, n,                          &
        &                                                n, n, n, n, y], [5, 5])
        integer, parameter :: sel(4) = [1, 3, 5, 7] !< S1 and the three components of T1.
        real(dp) :: h(2, 3), h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp), allocatable :: u(:, :), p(:, :), v(:, :), p_sf(:, :), w(:, :), e(:), e_sf(:)
        complex(dp), allocatable :: c(:), d(:)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: i

        h = reshape([1.0e-4_dp, 2.0e-4_dp, -3.0e-4_dp, 0.5e-4_dp, 0.0_dp, 4.0e-4_dp], [2, 3])
        call set_st_model(pot, 'components', 0.03_dp, h)
        call set_st_model(pot0, 'components', 0.03_dp, 0.0_dp * h)
        call set_st_model(pot_eff, 'effective_rss', 0.03_dp, h)
        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(sub, h0, h1, x0, x1)
        call sub%select_states([1, 1], [1, 3])
        call sub%update_geometry(reshape([0.1_dp, 0.2_dp], [1, 2]), 'model')
        call sub%eval()
        call setup_external(ext, 'true')

        ! Available transformations.
        ret = 1
        if (.not. pot%has_diabatic() .or. ext%has_diabatic()) return
        if (.not. transform_table_ok(pot, diab_table, .true.)) return
        if (.not. transform_table_ok(sub, diab_table, .true.)) return
        if (.not. transform_table_ok(pot_eff, diab_table, .true.)) return
        if (.not. transform_table_ok(ext, ext_table, .false.)) return
        if (pot%can_transform('mch', 'spin_free') .or. pot%can_transform('x', 'x')) return

        ! Composition, all states.
        ret = 2
        u = pot%get_transformation('diabatic', 'adiabatic')
        p = pot%get_transformation('diabatic', 'mch')
        v = pot%get_transformation('mch', 'adiabatic')
        p_sf = pot%get_transformation('diabatic_sf', 'mch_sf')
        if (any(u /= pot%adiab_trans) .or. any(u /= matmul(p, v))) return
        if (any(v /= pot%adiab_vec) .or. any(p_sf /= pot%mch_trans)) return
        if (any(p /= pot%model_states%expand(pot%mch_trans))) return
        if (any(abs(matmul(transpose(v), v) - unit_mat(7)) > tol)) return
        if (any(abs(matmul(transpose(p), p) - unit_mat(7)) > tol)) return
        if (any(pot%get_transformation('adiabatic', 'mch') /= transpose(v))) return
        if (any(pot%get_transformation('mch', 'diabatic') /= transpose(p))) return
        if (any(pot%get_transformation('mch_sf', 'diabatic_sf') /= transpose(p_sf))) return

        ! Composition, subset: S1 and T1 of two singlets and two triplets.
        ret = 3
        u = sub%get_transformation('diabatic', 'adiabatic')
        p = sub%get_transformation('diabatic', 'mch')
        v = sub%get_transformation('mch', 'adiabatic')
        p_sf = sub%get_transformation('diabatic_sf', 'mch_sf')
        if (any(shape(u) /= [8, 4]) .or. any(shape(p) /= [8, 4])) return
        if (any(shape(v) /= [4, 4]) .or. any(shape(p_sf) /= [4, 2])) return
        if (any(u /= sub%adiab_trans) .or. any(u /= matmul(p, v))) return
        if (any(p_sf /= sub%mch_trans(:, [1, 3]))) return
        w = sub%model_states%expand(sub%mch_trans)
        if (any(p /= w(:, sel))) return
        if (any(abs(matmul(transpose(v), v) - unit_mat(4)) > tol)) return
        if (any(abs(matmul(transpose(p), p) - unit_mat(4)) > tol)) return
        if (any(abs(matmul(transpose(p_sf), p_sf) - unit_mat(2)) > tol)) return

        ! Coefficients.
        ret = 4
        c = test_coefficients(7)
        u = pot%get_transformation('diabatic', 'adiabatic')
        v = pot%get_transformation('mch', 'adiabatic')
        d = pot%transform('adiabatic', 'diabatic', c)
        if (any(abs(d - matmul(u, c)) > tol)) return
        if (any(abs(pot%transform('diabatic', 'adiabatic', d) - c) > tol)) return
        d = pot%transform('adiabatic', 'mch', c)
        if (any(abs(d - matmul(v, c)) > tol)) return
        if (any(abs(pot%transform('mch', 'adiabatic', d) - c) > tol)) return
        c = test_coefficients(4)
        u = sub%get_transformation('diabatic', 'adiabatic')
        d = sub%transform('adiabatic', 'diabatic', c)
        if (size(d) /= 8 .or. any(abs(d - matmul(u, c)) > tol)) return
        if (any(abs(sub%transform('diabatic', 'adiabatic', d) - c) > tol)) return
        d = test_coefficients(8)
        c = sub%transform('diabatic', 'adiabatic', d)
        if (size(c) /= 4 .or. any(abs(c - matmul(transpose(u), d)) > tol)) return
        c = test_coefficients(4)
        d = sub%transform('diabatic_sf', 'mch_sf', c)
        if (size(d) /= 2) return

        ! Operators.
        ret = 5
        e = pot%get_energy('adiabatic')
        w = pot%transform('diabatic', 'adiabatic', pot%diab_w)
        if (any(shape(w) /= [7, 7]) .or. any(abs(w - diag(e)) > tol)) return
        w = pot%get_Hamiltonian('mch')
        e_sf = pot%get_energy('mch_sf')
        if (any(shape(w) /= [7, 7])) return
        if (any(abs(diag(w) - e_sf([1, 2, 3, 2, 3, 2, 3])) > tol)) return
        if (any(abs(w - pot%transform('diabatic', 'mch', pot%diab_w)) > tol)) return
        if (any(abs(pot%transform('mch', 'adiabatic', w) - diag(e)) > tol)) return
        w = pot%transform('diabatic_sf', 'mch_sf', pot%sf_w)
        if (any(abs(w - diag(e_sf)) > tol)) return
        if (any(abs(pot%get_Hamiltonian('mch_sf') - diag(e_sf)) > tol)) return
        if (any(pot%get_Hamiltonian('diabatic_sf') /= pot%sf_w)) return
        e = sub%get_energy('adiabatic')
        w = sub%transform('diabatic', 'adiabatic', sub%diab_w)
        if (any(shape(w) /= [4, 4]) .or. any(abs(w - diag(e)) > tol)) return
        if (any(abs(sub%transform('adiabatic', 'diabatic', w)                                   &
        &           - matmul(matmul(u, w), transpose(u))) > tol)) return
        w = sub%get_Hamiltonian('mch')
        if (any(shape(w) /= [4, 4])) return
        if (any(abs(sub%transform('mch', 'adiabatic', w) - diag(e)) > tol)) return
        if (any(abs(w - sub%transform('diabatic', 'mch', sub%diab_w)) > tol)) return

        ! Without spin-orbit couplings, each adiabatic state is one MCH state, with its energy
        ! (the adiabatic states are in the order of their energies, the MCH states in that of
        ! the spin components).
        ret = 6
        v = pot0%get_transformation('mch', 'adiabatic')
        e = pot0%get_energy('adiabatic')
        e_sf = pot0%get_energy('mch')
        do i = 1, 7
            if (count(abs(abs(v(:, i)) - 1.0_dp) < tol) /= 1) return
            if (count(abs(v(:, i)) < tol) /= 6) return
            if (abs(e_sf(maxloc(abs(v(:, i)), dim=1)) - e(i)) > tol) return
        end do
        if (any(abs(pot0%get_Hamiltonian('mch') - diag(pot0%get_energy('mch'))) > tol)) return

        ! Equal bases, and effective spin states.
        ret = 7
        c = test_coefficients(7)
        if (any(pot%transform('adiabatic', 'adiabatic', c) /= c)) return
        if (any(pot%transform('diabatic', 'diabatic', pot%diab_w) /= pot%diab_w)) return
        c = test_coefficients(4)
        if (any(sub%transform('mch', 'mch', c) /= c)) return
        c = test_coefficients(2)
        if (any(sub%transform('mch_sf', 'mch_sf', c) /= c)) return
        if (any(pot_eff%get_transformation('diabatic', 'mch')                                    &
        &       /= pot_eff%get_transformation('diabatic_sf', 'mch_sf'))) return
        ret = 0
    end function transform_test


    !> can_transform(bases(i), bases(j)) of p is table(i, j) for the bases diabatic, mch,
    !! adiabatic, diabatic_sf and mch_sf. With evaluated, get_transformation gives a matrix with
    !! the sizes of the two bases for each pair where it is true.
    logical function transform_table_ok(p, table, evaluated) result(ok)
        use evaluator_base_mod, only : potential_evaluator
        class(potential_evaluator), intent(in) :: p
        logical, intent(in) :: table(5, 5)
        logical, intent(in) :: evaluated
        character(len=11), parameter :: bases(5) = [character(len=11) :: 'diabatic', 'mch',      &
        &   'adiabatic', 'diabatic_sf', 'mch_sf']
        character(len=:), allocatable :: from, to
        integer :: i, j

        ok = .false.
        do j = 1, 5
            do i = 1, 5
                from = trim(bases(i))
                to = trim(bases(j))
                if (p%can_transform(from, to) .neqv. table(i, j)) return
                if (.not. table(i, j) .or. .not. evaluated) cycle
                if (any(shape(p%get_transformation(from, to))                                    &
                &       /= [p%basis_size(from), p%basis_size(to)])) return
            end do
        end do
        ok = .true.
    end function transform_table_ok


    !> Names of the bases (basis_mod) and the number of states of each basis: one singlet and two
    !! triplets with components and as effective states, S1 and T1 of two singlets and two
    !! triplets, and three singlets of an external program.
    integer function basis_test() result(ret)
        use basis_mod, only : basis_known, basis_is_spin_free, basis_full, basis_like
        use pot_test_mod, only : group_evaluator, linear_evaluator
        use external_json_mod, only : external_json_evaluator
        use global_defs, only : dp
        type(group_evaluator) :: pot, pot_eff
        type(linear_evaluator) :: sub
        type(external_json_evaluator) :: ext
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)

        ret = 1
        if (.not. (basis_known('diabatic') .and. basis_known('mch') .and. basis_known('adiabatic') &
        &          .and. basis_known('diabatic_sf') .and. basis_known('mch_sf'))) return
        if (basis_known('spin_free') .or. basis_known('adiabatic_sf') .or. basis_known('')) return
        if (.not. (basis_is_spin_free('diabatic_sf') .and. basis_is_spin_free('mch_sf'))) return
        if (basis_is_spin_free('diabatic') .or. basis_is_spin_free('mch')) return
        if (basis_is_spin_free('adiabatic')) return
        ret = 2
        if (basis_full('mch_sf') /= 'mch' .or. basis_full('diabatic_sf') /= 'diabatic') return
        if (basis_full('mch') /= 'mch' .or. basis_full('adiabatic') /= 'adiabatic') return
        if (basis_like('diabatic', 'mch_sf') /= 'diabatic_sf') return
        if (basis_like('diabatic', 'adiabatic') /= 'diabatic') return
        if (basis_like('diabatic_sf', 'mch') /= 'diabatic') return
        if (basis_like('mch', 'diabatic_sf') /= 'mch_sf') return
        if (basis_like('adiabatic', 'mch') /= 'adiabatic') return

        ret = 3
        call pot%set_model_states([1, 2], [1, 3], 'components')
        if (any(basis_sizes(pot) /= [7, 7, 7, 3, 3])) return
        call pot_eff%set_model_states([1, 2], [1, 3], 'effective_rss')
        if (any(basis_sizes(pot_eff) /= [3, 3, 3, 3, 3])) return
        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(sub, h0, h1, x0, x1)
        call sub%select_states([1, 1], [1, 3])
        if (any(basis_sizes(sub) /= [8, 4, 4, 4, 2])) return
        call setup_external(ext, 'true')
        if (any(basis_sizes(ext) /= [3, 3, 3, 3, 3])) return
        ret = 0
    end function basis_test


    !> Sizes of the bases diabatic, mch, adiabatic, diabatic_sf and mch_sf.
    function basis_sizes(p) result(n)
        use evaluator_base_mod, only : potential_evaluator
        class(potential_evaluator), intent(in) :: p
        integer :: n(5)

        n = [p%basis_size('diabatic'), p%basis_size('mch'), p%basis_size('adiabatic'),           &
        &    p%basis_size('diabatic_sf'), p%basis_size('mch_sf')]
    end function basis_sizes


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: transform_error
    !> @brief Transformations which should stop with an error.
    !> @details
    !! - transform_spaces - Between a basis with one state for each multiplet and a basis with
    !!   the spin components.
    !! - transform_no_diabatic - To the diabatic basis of an external program.
    !! - transform_wrong_point - From the states of a point which is not the previous point.
    !! - transform_no_link - From the states of another point, at a first point.
    !----------------------------------------------------------------------------------------------
    subroutine transform_error(case)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator, linear_evaluator
        use external_json_mod, only : external_json_evaluator
        character(len=*), intent(in) :: case
        type(group_evaluator) :: pot
        type(linear_evaluator) :: pot1, pot2, pot3
        type(external_json_evaluator) :: ext
        real(dp) :: h(2, 3), h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp), allocatable :: t(:, :)
        complex(dp), allocatable :: c(:)

        select case(case)
        case('transform_spaces')
            h = 1.0e-4_dp
            call set_st_model(pot, 'components', 0.03_dp, h)
            t = pot%get_transformation('mch_sf', 'adiabatic')
        case('transform_no_diabatic')
            call setup_external(ext, 'true')
            c = ext%transform('adiabatic', 'diabatic', test_coefficients(3))
        case('transform_wrong_point', 'transform_no_link')
            call general_couplings(h0, h1, x0, x1)
            call set_linear_model(pot1, h0, h1, x0, x1)
            call pot1%update_geometry(reshape([0.1_dp, 0.2_dp], [1, 2]), 'model')
            call pot1%eval()
            if (case == 'transform_no_link') then
                t = pot1%get_transformation('adiabatic', 'adiabatic', from_point=pot1%id + 100)
                return
            end if
            pot2 = pot1
            call pot2%eval(pot1)
            pot3 = pot2
            call pot3%eval(pot2)
            t = pot3%get_transformation('adiabatic', 'adiabatic', from_point=pot1%id)
        case default
            write(*, *) 'Unknown case: ', case
        end select
    end subroutine transform_error


    !> Normalised complex coefficients with all elements different from zero.
    function test_coefficients(n) result(c)
        use global_defs, only : dp
        integer, intent(in) :: n
        complex(dp) :: c(n)
        integer :: k

        do k = 1, n
            c(k) = cmplx(cos(0.7_dp * k) + 1.5_dp, sin(1.3_dp * k), kind=dp)
        end do
        c = c / norm2(abs(c))
    end function test_coefficients


    !> A point saved with save_state and loaded into a new evaluator gives the same state, and the
    !! next point evaluated with the loaded point as the previous one is the same as with the
    !! original (phases, link and, for effective spin states, the signs by continuity).
    !! - granucci-sst with all spin components.
    !! - Signed rss effective states of granucci-st, whose couplings change sign between the
    !!   points, so the sign references of the saved point decide the signs at the next point.
    !! - granucci-sst with the dynamics on one singlet and the triplet.
    integer function save_state_test() result(ret)
        use global_defs, only : dp
        use string_mod, only : string
        use model_1d_mod, only : model_1d_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(model_1d_evaluator) :: sst, st, sub
        type(effective_spin_evaluator) :: eff

        ret = 1
        call sst%init([string('granucci-sst')])
        if (.not. restart_same(sst, [9.0_dp, 9.5_dp, 10.0_dp])) return
        ret = 2
        call st%init([string('granucci-st'), string('r_sigma'), string('9.5')])
        call eff%init(st, 'rss', signed=.true.)
        if (.not. restart_same(eff, [9.0_dp, 10.0_dp, 10.5_dp])) return
        ret = 3
        call sub%init([string('granucci-sst')])
        call sub%select_states([1, 1], [1, 3])
        if (.not. restart_same(sub, [9.0_dp, 9.5_dp, 10.0_dp])) return
        ret = 0
    end function save_state_test


    !> Evaluate copies of template at x(1) and x(2), save the second point and load it into a new
    !! copy. Compare the saved states, and the states at x(3) evaluated from the original and
    !! from the loaded point, with the same id.
    logical function restart_same(template, x) result(same)
        use global_defs, only : dp
        use json_module, only : json_core, json_value
        use state_overlap_mod, only : last_requested_id, restore_last_id
        use evaluator_base_mod, only : potential_evaluator
        class(potential_evaluator), intent(in) :: template
        real(dp), intent(in) :: x(3)
        class(potential_evaluator), allocatable :: p1, p2, loaded, n1, n2
        type(json_core) :: json
        type(json_value), pointer :: root
        character(len=:), allocatable :: text_saved, text_loaded, text_n1, text_n2
        integer :: id0

        same = .false.
        allocate(p1, source=template)
        allocate(p2, source=template)
        allocate(loaded, source=template)
        allocate(n1, source=template)
        allocate(n2, source=template)
        call p1%update_geometry(reshape([x(1)], [1, 1]), 'model')
        call p1%eval()
        call p2%update_geometry(reshape([x(2)], [1, 1]), 'model')
        call p2%eval(p1)
        text_saved = state_text(p2)

        call json%initialize()
        call json%deserialize(root, text_saved)
        call loaded%load_state(json, root)
        call json%destroy(root)
        text_loaded = state_text(loaded)
        if (text_loaded /= text_saved) return

        id0 = last_requested_id()
        call n1%update_geometry(reshape([x(3)], [1, 1]), 'model')
        call n1%eval(p2)
        text_n1 = state_text(n1)
        call restore_last_id(id0)
        call n2%update_geometry(reshape([x(3)], [1, 1]), 'model')
        call n2%eval(loaded)
        text_n2 = state_text(n2)
        same = (text_n1 == text_n2)
    end function restart_same


    !> The state of a point, as JSON text.
    function state_text(pot) result(text)
        use json_module, only : json_core, json_value
        use evaluator_base_mod, only : potential_evaluator
        class(potential_evaluator), intent(in) :: pot
        character(len=:), allocatable :: text
        type(json_core) :: json
        type(json_value), pointer :: root

        call json%initialize()
        call json%create_object(root, '')
        call pot%save_state(json, root)
        call json%serialize(root, text)
        call json%destroy(root)
    end function state_text


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
        call pot%set_model_states([1, 2], [1, 3], 'components')
        if (any(pot%state_map /= [1, 2, 5, 3, 6, 4, 7])) return
        call pot_zh%set_model_states([1, 2], [1, 3], 'components')
        if (any(pot_zh%state_map /= [(i, i = 1, 7)])) return
        pot_eff%component_order = 'by_state'
        call pot_eff%set_model_states([1, 2], [1, 3], 'effective_rss')
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
    !! flat, and the spin-orbit couplings within groups are checked at the first evaluation.
    integer function derivative_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : coord_evaluator
        type(coord_evaluator) :: pot
        real(dp), allocatable :: grad(:), nadv(:), d(:)
        real(dp) :: geom(3, 2)

        call pot%set_model_states([2], [1], 'components')
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
        ret = 0
    end function derivative_test


    !----------------------------------------------------------------------------------------------
    !> The link between the states of two points. Each evaluation gets a new id. With a previous
    !! point, the phases of the states are matched to it and the link holds U(prev)^T U and its
    !! orthogonalisation, which is the same for exact eigenvectors. Without a previous point, a
    !! copied link is removed. A point evaluated from another previous point does not link the
    !! first one.
    !----------------------------------------------------------------------------------------------
    integer function link_test() result(ret)
        use global_defs, only : dp
        use matrix_mod, only : unit_mat
        use evaluator_base_mod, only : new_request
        use pot_test_mod, only : linear_evaluator
        type(linear_evaluator) :: pot1, pot2, pot3
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2), x(1, 2)
        real(dp), allocatable :: s(:, :), t(:, :), u(:, :)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: i

        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(pot1, h0, h1, x0, x1)
        call pot1%update_geometry(reshape([0.1_dp, 0.2_dp], [1, 2]), 'model')
        call pot1%eval()
        ret = 1
        if (pot1%id <= 0 .or. pot1%link%valid() .or. pot1%previous_point() /= 0) return
        if (pot1%id_diab <= 0 .or. pot1%id_diab == pot1%id) return
        if (.not. pot1%diab%matches(pot1%id_diab, 'diabatic_sf', pot1%id, 'mch_sf')) return

        ! The next point, as in the trajectory data a copy of the previous one, after a small step.
        ret = 2
        pot2 = pot1
        call pot2%update_geometry(reshape([0.1_dp + 1.0e-6_dp, 0.2_dp - 1.0e-6_dp], [1, 2]),      &
        &                         'model')
        call pot2%eval(pot1)
        if (pot2%id == pot1%id) return
        if (.not. pot2%link%matches(pot1%id, 'adiabatic', pot2%id, 'adiabatic')) return
        if (pot2%link%matches(pot1%id, 'mch', pot2%id, 'mch')) return
        if (pot2%id_diab /= pot1%id_diab) return
        if (.not. pot2%diab%matches(pot1%id_diab, 'diabatic_sf', pot2%id, 'mch_sf')) return
        s = matmul(transpose(pot1%adiab_trans), pot2%adiab_trans)
        if (any(abs(pot2%link%s - s) > tol)) return
        if (any(abs(pot2%link%t - s) > tol)) return
        if (any(abs(matmul(transpose(pot2%link%t), pot2%link%t) - unit_mat(8)) > tol)) return
        ! After a small step, each state overlaps mostly with itself, with a matched phase.
        do i = 1, 8
            if (pot2%link%s(i, i) < 0.99_dp) return
        end do

        ! Without a previous point, the copied link is removed.
        ret = 3
        pot3 = pot2
        call pot3%eval()
        if (pot3%link%valid() .or. pot3%id == pot2%id) return

        ! The next point links the point it was evaluated from, not the one before it.
        ret = 4
        pot3 = pot2
        call pot3%update_geometry(reshape([0.1_dp + 2.0e-6_dp, 0.2_dp], [1, 2]), 'model')
        call pot3%eval(pot2)
        if (.not. pot3%link%matches(pot2%id, 'adiabatic', pot3%id, 'adiabatic')) return
        if (pot3%link%matches(pot1%id, 'adiabatic', pot3%id, 'adiabatic')) return

        ! The states of the previous point through get_transformation: the orthogonalised
        ! overlaps, their transpose in the other direction, and followed by the transformation
        ! to another basis of this point. The id of the point itself gives its own states.
        ret = 5
        if (pot2%previous_point() /= pot1%id) return
        t = pot2%get_transformation('adiabatic', 'adiabatic', from_point=pot1%id)
        if (any(t /= pot2%link%t)) return
        if (any(pot2%get_transformation('adiabatic', 'adiabatic', to_point=pot1%id)              &
        &       /= transpose(t))) return
        u = pot2%get_transformation('diabatic', 'adiabatic')
        if (any(abs(pot2%get_transformation('adiabatic', 'diabatic', from_point=pot1%id)         &
        &           - matmul(t, transpose(u))) > tol)) return
        if (any(abs(pot2%get_transformation('diabatic', 'adiabatic', to_point=pot1%id)           &
        &           - matmul(u, transpose(t))) > tol)) return
        if (any(pot2%get_transformation('diabatic', 'adiabatic', from_point=pot2%id,             &
        &                               to_point=pot2%id) /= u)) return
        if (any(abs(pot2%transform('adiabatic', 'adiabatic', pot1%adiab_w, from_point=pot1%id)   &
        &           - matmul(matmul(transpose(t), pot1%adiab_w), t)) > tol)) return

        ! The link is in the basis of the request (the basis of the dynamics): the overlaps of
        ! the MCH vectors of each manifold for mch_sf, the same for every component for mch, and
        ! the unit matrix for the diabatic states, which are the same at all points.
        ret = 6
        x = reshape([0.1_dp + 1.0e-6_dp, 0.2_dp - 1.0e-6_dp], [1, 2])
        pot3 = pot1
        call pot3%set_request(new_request('mch_sf'))
        call pot3%update_geometry(x, 'model')
        call pot3%eval(pot1)
        if (.not. pot3%link%matches(pot1%id, 'mch_sf', pot3%id, 'mch_sf')) return
        s = matmul(transpose(pot1%mch_trans), pot3%mch_trans)
        if (any(shape(pot3%link%s) /= [4, 4]) .or. any(abs(pot3%link%s - s) > tol)) return
        if (any(pot3%link%s(1:2, 3:4) /= 0.0_dp)) return
        if (any(pot3%get_transformation('mch_sf', 'mch_sf', from_point=pot1%id)                  &
        &       /= pot3%link%t)) return
        if (any(abs(pot3%adiab_w - pot2%adiab_w) > tol) .or.                                     &
        &   any(abs(pot3%adiab_trans - pot2%adiab_trans) > tol)) return
        ret = 7
        pot3 = pot1
        call pot3%set_request(new_request('mch'))
        call pot3%update_geometry(x, 'model')
        call pot3%eval(pot1)
        if (.not. pot3%link%matches(pot1%id, 'mch', pot3%id, 'mch')) return
        if (any(abs(pot3%link%s - pot1%model_states%expand(s)) > tol)) return
        ret = 8
        pot3 = pot1
        call pot3%set_request(new_request('diabatic'))
        call pot3%update_geometry(x, 'model')
        call pot3%eval(pot1)
        if (.not. pot3%link%matches(pot1%id, 'diabatic', pot3%id, 'diabatic')) return
        if (any(pot3%link%s /= unit_mat(8)) .or. any(pot3%link%t /= unit_mat(8))) return
        ret = 0
    end function link_test


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
        call pot%set_model_states([4, 3], [1, 3], 'components')
        if (pot%model_states%n_state /= 13 .or. pot%model_states%n_sf /= 7) return
        if (size(pot%model_states%manifold) /= 2) return
        if (any(pot%model_states%manifold%n /= [4, 3])) return
        if (any(pot%model_states%manifold%n_comp /= [1, 3])) return
        if (any(pot%model_states%manifold%mult /= [1, 3])) return
        if (any(pot%model_states%manifold%i0 /= [1, 5])) return
        if (any(pot%model_states%manifold%sf0 /= [1, 5])) return

        ret = 2
        loc = pot%model_states%locate_state(4)
        if (loc%group /= 1 .or. loc%manifold /= 1 .or. loc%mult /= 1) return
        if (loc%component /= 1 .or. loc%index /= 4 .or. loc%sf /= 4) return
        loc = pot%model_states%locate_state(5)
        if (loc%group /= 2 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 1 .or. loc%index /= 1 .or. loc%sf /= 5) return
        loc = pot%model_states%locate_state(10)
        if (loc%group /= 3 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 2 .or. loc%index /= 3 .or. loc%sf /= 7) return
        loc = pot%model_states%locate_state(13)
        if (loc%group /= 4 .or. loc%manifold /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 3 .or. loc%index /= 3 .or. loc%sf /= 7) return

        ! The mch_sf basis has one state for each spin-free state. Summed over the components,
        ! a matrix of ones counts the component pairs, and a vector of ones the components.
        ret = 3
        if (pot%basis_size('mch_sf') /= 7 .or. pot%basis_size('diabatic') /= 13) return
        sf_ones = pot%model_states%sum_components(ones)
        if (sf_ones(1, 4) /= 1.0_dp .or. sf_ones(2, 5) /= 3.0_dp .or. sf_ones(7, 3) /= 3.0_dp) return
        if (sf_ones(5, 7) /= 9.0_dp .or. sum(sf_ones) /= 13.0_dp**2) return
        if (pot%basis_size('mch') /= 13 .or. pot%basis_size('diabatic_sf') /= 7) return
        if (any(pot%model_states%sum_components(ones(:, 1))                                      &
        &       /= [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 3.0_dp])) return
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
        if (pot%model_states%n_state /= 8 .or. pot%model_states%n_sf /= 4) return
        if (any(pot%model_states%manifold%n /= [2, 2])) return
        if (any(pot%model_states%manifold%n_comp /= [1, 3])) return
        if (any(pot%model_states%manifold%i0 /= [1, 3])) return
        if (any(pot%model_states%manifold%mult /= [1, 3])) return
        loc = pot%model_states%locate_state(6)
        if (loc%manifold /= 2 .or. loc%component /= 2 .or. loc%index /= 2) return

        ret = 2
        if (pot%model_states%max_soc_in_groups(model%soc%zero) /= 0.0_dp) return

        ! Coupling between the two states of the second triplet component.
        ret = 3
        soc(5, 6) = 1.0e-4_dp
        soc(6, 5) = 1.0e-4_dp
        call model%soc%add_zero(model%nmode, soc)
        call pot_bad%initialize(model)
        if (pot_bad%model_states%max_soc_in_groups(model%soc%zero) <= 0.0_dp) return
        ret = 0
    end function vc_group_test


    !> The mch_sf and mch bases of one singlet and two triplets, with couplings h(k, c) between
    !! the singlet and component c of triplet k, and the coupling mix between the two triplets.
    integer function mch_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        type(group_evaluator) :: pot
        type(group_evaluator) :: pot_eff
        real(dp) :: h(2, 3), r(3, 3), s2(3, 3), s2_ref(3, 3), u(3, 3)
        real(dp), allocatable :: s2_c(:, :), w(:, :), e_c(:)
        real(dp), parameter :: tiny = 1.0e-14_dp
        integer, parameter :: sf(7) = [1, 2, 3, 2, 3, 2, 3] !< Spin-free state of each state.
        integer :: k

        h = reshape([1.0e-4_dp, 2.0e-4_dp, -3.0e-4_dp, 0.5e-4_dp, 0.0_dp, 4.0e-4_dp], [2, 3])

        ! Without mixing, the squared coupling of the singlet and triplet k is sum_c h(k, c)**2.
        ret = 1
        call set_st_model(pot, 'components', 0.0_dp, h)
        if (pot%basis_size('mch_sf') /= 3 .or. pot%basis_size('adiabatic') /= 7) return
        if (size(pot%get_energy('mch_sf')) /= 3) return
        if (any(abs(pot%get_energy('mch_sf') - [0.0_dp, 0.1_dp, 0.2_dp]) > tiny)) return
        s2 = pot%get_soc2('mch_sf')
        do k = 1, 2
            if (abs(s2(1, k + 1) - sum(h(k, :)**2)) > tiny) return
        end do
        if (any(abs(s2 - transpose(s2)) > tiny)) return

        ! With mixing of the triplets, the couplings are those of the spin-free states, and their
        ! sum over the triplets does not change.
        ret = 2
        call set_st_model(pot, 'components', 0.03_dp, h)
        s2 = pot%get_soc2('mch_sf')
        u = pot%get_transformation('diabatic_sf', 'mch_sf')
        if (any(shape(u) /= [3, 3]) .or. any(u /= pot%mch_trans)) return
        do k = 1, 2
            if (abs(s2(1, k + 1) - sum(matmul(u(2:3, k + 1), h)**2)) > tiny) return
        end do
        if (abs(s2(1, 2) + s2(1, 3) - sum(h**2)) > tiny) return
        if (abs(pot%get_energy('mch_sf', 3) - pot%get_energy('mch_sf', 2)                  &
        &       - 2 * sqrt(0.05_dp**2 + 0.03_dp**2)) > tiny) return

        ! A rotation of the triplet components does not change the couplings.
        ret = 3
        s2_ref = s2
        r = reshape([0.6_dp, 0.8_dp, 0.0_dp, -0.8_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp],    &
        &           [3, 3])
        r = matmul(r, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.28_dp, 0.96_dp, 0.0_dp,          &
        &                      -0.96_dp, 0.28_dp], [3, 3]))
        call set_st_model(pot, 'components', 0.03_dp, matmul(h, r))
        if (any(abs(pot%get_soc2('mch_sf') - s2_ref) > tiny)) return

        ! Effective states with the root-sum-square couplings give the same couplings without
        ! mixing, and with mixing if the component couplings of the two triplets are parallel.
        ret = 4
        call set_st_model(pot, 'components', 0.0_dp, h)
        call set_st_model(pot_eff, 'effective_rss', 0.0_dp, h)
        if (pot_eff%basis_size('mch_sf') /= 3 .or. pot_eff%model_states%n_state /= 3) return
        if (any(abs(pot_eff%get_soc2('mch_sf') - pot%get_soc2('mch_sf')) > tiny)) return
        h(2, :) = 0.5_dp * h(1, :)
        call set_st_model(pot, 'components', 0.03_dp, h)
        call set_st_model(pot_eff, 'effective_rss', 0.03_dp, h)
        if (any(abs(pot_eff%get_soc2('mch_sf') - pot%get_soc2('mch_sf')) > tiny)) return

        ! The mch basis has the spin components of the MCH states: the energies of the spin-free
        ! states, the spin-orbit couplings off the diagonal of its Hamiltonian, and their squares
        ! summed over the components for the mch_sf basis.
        ret = 5
        h = reshape([1.0e-4_dp, 2.0e-4_dp, -3.0e-4_dp, 0.5e-4_dp, 0.0_dp, 4.0e-4_dp], [2, 3])
        call set_st_model(pot, 'components', 0.03_dp, h)
        e_c = pot%get_energy('mch')
        if (size(e_c) /= 7) return
        do k = 1, 7
            if (e_c(k) /= pot%get_energy('mch_sf', sf(k))) return
            if (pot%get_energy('mch', k) /= e_c(k)) return
        end do
        w = pot%get_Hamiltonian('mch')
        s2_c = pot%get_soc2('mch')
        do k = 1, 7
            if (w(k, k) /= e_c(k)) return
            w(k, k) = 0.0_dp
        end do
        if (any(shape(s2_c) /= [7, 7]) .or. any(abs(s2_c - w**2) > tiny)) return
        if (any(abs(pot%states%sum_components(s2_c) - pot%get_soc2('mch_sf')) > tiny)) return
        ret = 0
    end function mch_test


    !> Set and evaluate the potential of mch_test. With effective states, the coupling of
    !! the singlet and triplet k is the root-sum-square of h(k, :).
    subroutine set_st_model(pot, spin_states, mix, h)
        use global_defs, only : dp
        use pot_test_mod, only : group_evaluator
        type(group_evaluator), intent(inout) :: pot
        character(len=*), intent(in) :: spin_states
        real(dp), intent(in) :: mix
        real(dp), intent(in) :: h(:, :)
        integer :: k, c, i

        if (.not. pot%model_states%is_set()) then
            call pot%set_model_states([1, 2], [1, 3], spin_states)
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
        associate(ms => eff%model_states)
            if (ms%n_state /= 4 .or. ms%n_sf /= 4 .or. ms%spin_states /= 'effective_rss') return
            if (any(ms%manifold%n_comp /= [1, 1]) .or. any(ms%manifold%mult /= [1, 3])) return
        end associate
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
        if (any(abs(eff%get_energy('mch_sf') - src%get_energy('mch_sf')) > tol)) return
        if (allocated(eff%source%adiab_w) .or. allocated(eff%source%mch_trans)) return
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
        e_sf = src%get_energy('mch_sf')
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
                call eff_d%eval(eff(m))
                e_p = eff_d%get_energy('adiabatic')
                geom_d(1, n) = geom(1, n) - delta
                call eff_d%update_geometry(geom_d, 'model')
                call eff_d%eval(eff(m))
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
    !! Each point is evaluated with the previous one. A second evaluation at the same geometry, with
    !! the same previous point or with the point itself as the previous one, gives the same
    !! couplings.
    !----------------------------------------------------------------------------------------------
    integer function effective_sign_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff(5)
        type(effective_spin_evaluator) :: prev !< Previous point.
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
                if (i == 0) then
                    call eff(m)%eval()
                else
                    call eff(m)%eval(prev)
                end if
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
                ! A second evaluation, with the same previous point and with the current point as
                ! the previous one.
                if (i == 5 .or. i == 6) then
                    soc_w = eff(m)%soc_w
                    call eff(m)%eval(prev)
                    if (any(eff(m)%soc_w /= soc_w)) return
                    prev = eff(m)
                    call eff(m)%eval(prev)
                    if (any(eff(m)%soc_w /= soc_w)) return
                end if
                prev = eff(m)
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
        if (st%model_states%n_state /= 4 .or. st%model_states%n_sf /= 2) return
        if (any(st%model_states%manifold%mult /= [1, 3])) return
        if (st%model_states%spin_states /= 'components') return

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
        if (sst%model_states%n_state /= 5 .or. sst%model_states%n_sf /= 3) return
        if (any(sst%model_states%manifold%n /= [2, 1])) return
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


    !> The state layout of a potential with components, of a potential with effective states,
    !! of the effective model built from the first, which adds the model and the layout of its
    !! source, and of a potential with the dynamics on a subset of its states.
    integer function state_layout_test() result(ret)
        use string_mod, only : string
        use json_module, only : json_core, json_value
        use pot_test_mod, only : group_evaluator
        use model_1d_mod, only : model_1d_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(model_1d_evaluator) :: st, sst
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
        if (jint(json, root, 'n_state.mch') /= 4) return
        if (jint(json, root, 'n_state.mch_sf') /= 2) return
        if (jint(json, root, 'n_state.diabatic_sf') /= 2) return
        if (jint(json, root, 'manifolds(1).multiplicity') /= 1) return
        if (jint(json, root, 'manifolds(1).n') /= 1) return
        if (jint(json, root, 'manifolds(2).multiplicity') /= 3) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 3) return
        if (jint(json, root, 'manifolds(2).first_spin_free') /= 2) return
        if (jint(json, root, 'manifolds(2).first_state') /= 2) return
        call json%destroy(root)

        ! Effective states given as input.
        ret = 2
        call grp%set_model_states([1, 2], [1, 3], 'effective_rss')
        call json%create_object(root, '')
        call grp%add_state_layout(json, root)
        if (jstr(json, root, 'spin_states') /= 'effective_rss') return
        if (jint(json, root, 'n_state.diabatic') /= 3) return
        if (jint(json, root, 'n_state.mch') /= 3) return
        if (jint(json, root, 'n_state.mch_sf') /= 3) return
        if (jint(json, root, 'manifolds(2).n') /= 2) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 1) return
        call json%destroy(root)

        ! Effective model: its own layout, the model, and the layout of the source.
        ret = 3
        call eff%init(st, 'sq3', signed=.false.)
        call json%create_object(root, '')
        call eff%add_state_layout(json, root)
        if (jstr(json, root, 'spin_states') /= 'effective_rss') return
        if (jint(json, root, 'n_state.diabatic') /= 2) return
        if (jint(json, root, 'manifolds(2).n_comp') /= 1) return
        if (jstr(json, root, 'effective_spin_states.model') /= 'sq3') return
        if (jlog(json, root, 'effective_spin_states.signed')) return
        if (jstr(json, root, 'source.spin_states') /= 'components') return
        if (jint(json, root, 'source.n_state.diabatic') /= 4) return
        if (jint(json, root, 'source.manifolds(2).n_comp') /= 3) return
        call json%destroy(root)

        ! Subset of the states: the states of the dynamics, and the states of the model.
        ret = 4
        call sst%init([string('granucci-sst')])
        call sst%select_states([1, 1], [1, 3])
        call json%create_object(root, '')
        call sst%add_state_layout(json, root)
        if (jint(json, root, 'n_state.diabatic') /= 5) return
        if (jint(json, root, 'n_state.adiabatic') /= 4) return
        if (jint(json, root, 'n_state.mch') /= 4) return
        if (jint(json, root, 'n_state.mch_sf') /= 2) return
        if (jint(json, root, 'n_state.diabatic_sf') /= 3) return
        if (jint(json, root, 'manifolds(1).n') /= 1) return
        if (jint(json, root, 'manifolds(2).first_state') /= 2) return
        if (jstr(json, root, 'model.spin_states') /= 'components') return
        if (jint(json, root, 'model.manifolds(1).n') /= 2) return
        if (jint(json, root, 'model.manifolds(2).first_state') /= 3) return
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


    !> Subsets of a state_layout: the first spin-free states of each multiplicity, in every
    !! component, and their positions in the full layout. A matrix with the rows over the
    !! spin-free states of the full layout and the columns over those of the subset is expanded
    !! to the selected columns of the expanded square matrix.
    integer function layout_subset_test() result(ret)
        use global_defs, only : dp
        use state_layout_mod, only : state_layout
        type(state_layout) :: full, sub
        integer, allocatable :: idx(:), idx_sf(:), rest(:)
        real(dp) :: sq(7, 7)
        real(dp), allocatable :: m(:, :), m_sub(:, :)
        integer :: i, j

        ! Four singlets and three triplets: S1-S4, then T1-T3 for each component.
        ret = 1
        call full%set([4, 3], [1, 3], 'components')
        sub = full%subset([2, 1], [1, 3])
        if (sub%n_state /= 5 .or. sub%n_sf /= 3 .or. sub%spin_states /= 'components') return
        if (any(sub%manifold%n /= [2, 1]) .or. any(sub%manifold%i0 /= [1, 3])) return
        call full%subset_index(sub, idx, idx_sf, rest)
        if (any(idx /= [1, 2, 5, 8, 11]) .or. any(idx_sf /= [1, 2, 5])) return
        if (any(rest /= [3, 4, 6, 7, 9, 10, 12, 13])) return
        do j = 1, 7
            do i = 1, 7
                sq(i, j) = real(i + 10 * j, dp)
            end do
        end do
        m = full%expand(sq)
        m_sub = full%expand(sq(:, idx_sf), sub)
        if (any(shape(m) /= [13, 13]) .or. any(shape(m_sub) /= [13, 5])) return
        if (any(m_sub /= m(:, idx))) return
        if (any(sub%expand(sq(idx_sf, :), full) /= m(idx, :))) return

        ! The manifolds keep the order of the full layout.
        ret = 2
        sub = full%subset([1, 2], [3, 1])
        if (any(sub%manifold%mult /= [1, 3]) .or. any(sub%manifold%n /= [2, 1])) return

        ! A multiplicity which is not given is dropped.
        ret = 3
        sub = full%subset([2], [3])
        if (size(sub%manifold) /= 1 .or. sub%n_state /= 6 .or. sub%n_sf /= 2) return
        call full%subset_index(sub, idx, idx_sf, rest)
        if (any(idx /= [5, 6, 8, 9, 11, 12]) .or. any(idx_sf /= [5, 6])) return
        if (size(rest) /= 7) return
        m_sub = full%expand(sq(:, idx_sf), sub)
        if (any(shape(m_sub) /= [13, 6]) .or. any(m_sub /= m(:, idx))) return

        ! All states, and effective states.
        ret = 4
        sub = full%subset([4, 3], [1, 3])
        call full%subset_index(sub, idx, idx_sf, rest)
        if (any(idx /= [(i, i = 1, 13)]) .or. size(rest) /= 0) return
        call full%set([4, 3], [1, 3], 'effective_rss')
        sub = full%subset([2, 1], [1, 3])
        call full%subset_index(sub, idx, idx_sf, rest)
        if (any(idx /= [1, 2, 5]) .or. any(idx_sf /= [1, 2, 5])) return
        if (any(rest /= [3, 4, 6, 7])) return
        ret = 0
    end function layout_subset_test


    !> Dynamics on a subset of the states of the linear model (S1 and T1 of two singlets and two
    !! triplets): sizes of the bases, the adiabatic states as eigenstates of the Hamiltonian in
    !! the selected MCH states, and the mch basis as a part of that of the full model.
    integer function subset_test() result(ret)
        use global_defs, only : dp
        use linalg_wrapper_mod, only : syev
        use matrix_mod, only : unit_mat
        use pot_test_mod, only : linear_evaluator
        type(linear_evaluator) :: full, every, sub, nosoc
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: x(1, 2), e(4)
        real(dp), allocatable :: p(:, :), h_sub(:, :), u(:, :), e_mch(:), s2(:, :), tr(:, :)
        integer, parameter :: sel(4) = [1, 3, 5, 7] !< S1 and the three components of T1.
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: n

        x = reshape([0.1_dp, 0.2_dp], [1, 2])
        call general_couplings(h0, h1, x0, x1)
        call set_linear_model(full, h0, h1, x0, x1)
        call full%update_geometry(x, 'model')
        call full%eval()

        ! All states selected: the same as without a selection.
        ret = 1
        call set_linear_model(every, h0, h1, x0, x1)
        call every%select_states([2, 2], [3, 1])
        call every%update_geometry(x, 'model')
        call every%eval()
        if (any(every%adiab_trans /= full%adiab_trans)) return
        if (any(every%adiab_w /= full%adiab_w) .or. any(every%adiab_dw /= full%adiab_dw)) return

        ! Sizes of the bases with S1 and T1 selected.
        ret = 2
        call set_linear_model(sub, h0, h1, x0, x1)
        call sub%select_states([1, 1], [1, 3])
        call sub%update_geometry(x, 'model')
        call sub%eval()
        if (sub%basis_size('diabatic') /= 8 .or. sub%basis_size('adiabatic') /= 4) return
        if (sub%basis_size('mch_sf') /= 2) return
        if (any(sub%sel /= sel) .or. any(sub%sel_mch /= [1, 3])) return
        if (any(sub%uns /= [2, 4, 6, 8])) return
        u = sub%get_transformation('diabatic', 'adiabatic')
        if (any(shape(u) /= [8, 4])) return
        if (any(shape(sub%get_transformation('adiabatic', 'diabatic')) /= [4, 8])) return
        if (any(shape(sub%get_transformation('adiabatic', 'adiabatic')) /= [4, 4])) return
        if (any(shape(sub%get_transformation('diabatic', 'diabatic')) /= [8, 8])) return
        if (any(shape(sub%get_Hamiltonian('adiabatic')) /= [4, 4])) return
        if (any(shape(sub%get_Hamiltonian('diabatic')) /= [8, 8])) return
        if (size(sub%get_energy('adiabatic')) /= 4) return
        if (size(sub%get_energy('diabatic')) /= 8) return
        if (size(sub%get_gradient('adiabatic', 4)) /= 2) return

        ! The adiabatic states are orthonormal, lie in the selected MCH states of the model, and
        ! their energies are the eigenvalues of the Hamiltonian in these states. They are not
        ! eigenstates of the full model.
        ret = 3
        if (any(abs(matmul(transpose(u), u) - unit_mat(4)) > tol)) return
        p = full%model_states%expand(full%mch_trans)
        if (any(abs(matmul(p(:, sel), matmul(transpose(p(:, sel)), u)) - u) > tol)) return
        h_sub = matmul(transpose(p(:, sel)), matmul(full%diab_w, p(:, sel)))
        call syev(h_sub, e, jobz='N', uplo='U')
        if (any(abs(sub%get_energy('adiabatic') - e) > tol)) return
        if (all(abs(e - full%adiab_w(1, 1)) < 1.0e-9_dp)) return
        if (minval(abs(e - full%adiab_w(1, 1))) < 1.0e-9_dp) return

        ! The mch basis is a part of that of the full model.
        ret = 4
        e_mch = full%get_energy('mch_sf')
        if (any(sub%get_energy('mch_sf') /= e_mch([1, 3]))) return
        if (sub%get_energy('mch_sf', 2) /= e_mch(3)) return
        if (any(sub%get_gradient('mch_sf', 2) /= full%get_gradient('mch_sf', 3))) return
        if (any(sub%get_nadv('mch_sf', 1, 2) /= full%get_nadv('mch_sf', 1, 3))) return
        s2 = full%get_soc2('mch_sf')
        if (any(sub%get_soc2('mch_sf') /= s2([1, 3], [1, 3]))) return
        tr = full%get_transformation('diabatic_sf', 'mch_sf')
        if (any(shape(sub%get_transformation('diabatic_sf', 'mch_sf')) /= [4, 2])) return
        if (any(sub%get_transformation('diabatic_sf', 'mch_sf') /= tr(:, [1, 3]))) return
        if (any(shape(sub%get_soc2('diabatic')) /= [8, 8])) return

        ! The mch basis: S1 and the three components of T1 (full model: S1, S2, then T1 and T2
        ! for each component). Gradients and couplings are those of the spin-free states, with
        ! no coupling between different components or multiplicities.
        ret = 5
        if (sub%basis_size('mch') /= 4 .or. full%basis_size('mch') /= 8) return
        if (any(sub%get_energy('mch') /= e_mch([1, 3, 3, 3]))) return
        if (any(full%get_energy('mch') /= e_mch([1, 2, 3, 4, 3, 4, 3, 4]))) return
        if (any(sub%get_gradient('mch', 3) /= full%get_gradient('mch_sf', 3))) return
        if (any(full%get_gradient('mch', 6) /= full%get_gradient('mch_sf', 4))) return
        if (any(full%get_nadv('mch', 5, 6) /= full%get_nadv('mch_sf', 3, 4))) return
        if (any(full%get_nadv('mch', 1, 2) /= full%get_nadv('mch_sf', 1, 2))) return
        if (all(full%get_nadv('mch', 5, 6) == 0.0_dp)) return
        if (any(full%get_nadv('mch', 3, 6) /= 0.0_dp)) return
        if (any(full%get_nadv('mch', 2, 3) /= 0.0_dp)) return
        if (any(full%get_nadv('diabatic_sf', 1, 2) /= 0.0_dp)) return
        if (any(full%get_gradient('diabatic_sf', 3) /= full%sf_dw(:, 3, 3))) return
        if (any(full%get_energy('diabatic_sf') /= [(full%sf_w(n, n), n = 1, 4)])) return
        s2 = full%get_soc2('mch')
        if (any(sub%get_soc2('mch') /= s2(sel, sel))) return
        if (any(abs(sub%states%sum_components(sub%get_soc2('mch')) - sub%get_soc2('mch_sf'))     &
        &       > tol)) return

        ! Without spin-orbit couplings, the adiabatic states are the selected MCH states: T1 three
        ! times and S1, with the gradients of the MCH states.
        ret = 6
        h0 = 0.0_dp
        h1 = 0.0_dp
        x0 = 0.0_dp
        x1 = 0.0_dp
        call set_linear_model(nosoc, h0, h1, x0, x1)
        call nosoc%select_states([1, 1], [1, 3])
        call nosoc%update_geometry(x, 'model')
        call nosoc%eval()
        e_mch = nosoc%get_energy('mch_sf')
        if (e_mch(1) >= e_mch(2)) return
        if (any(abs(nosoc%get_energy('adiabatic') - e_mch([1, 2, 2, 2])) > tol)) return
        if (any(abs(nosoc%get_gradient('adiabatic', 1) - nosoc%get_gradient('mch_sf', 1))        &
        &       > tol)) return
        do n = 2, 4
            if (any(abs(nosoc%get_gradient('adiabatic', n) - nosoc%get_gradient('mch_sf', 2))    &
            &       > tol)) return
        end do
        ret = 0
    end function subset_test


    !> Gradients and nonadiabatic couplings of the adiabatic states of a subset (S1 and T1 of the
    !! linear model), against central differences of the energies and of the overlaps with the
    !! states of displaced points. With exact_gradient they agree. Without it, the term from the
    !! change of the selected MCH states is missing, and the gradients differ.
    integer function subset_gradient_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        type(linear_evaluator) :: pot, pot_hf, pp, pm
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: x(1, 2), xd(1, 2), fd(4, 2), nac_fd(4, 4, 2), err, err_hf
        real(dp), allocatable :: g(:), nac(:)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp), parameter :: tol = 1.0e-8_dp
        integer :: n, a, b

        x = reshape([0.1_dp, 0.2_dp], [1, 2])
        call general_couplings(h0, h1, x0, x1)
        ! Larger couplings, so that the states of T1 are well separated and the missing term is
        ! far above the error of the central differences.
        h0 = 5.0_dp * h0
        h1 = 5.0_dp * h1
        x0 = 5.0_dp * x0
        x1 = 5.0_dp * x1
        call set_linear_model(pot, h0, h1, x0, x1)
        call pot%select_states([1, 1], [1, 3])
        call pot%update_geometry(x, 'model')
        call pot%eval()

        do n = 1, 2
            pp = pot
            xd = x
            xd(1, n) = x(1, n) + step
            call pp%update_geometry(xd, 'model')
            call pp%eval(pot)
            pm = pot
            xd(1, n) = x(1, n) - step
            call pm%update_geometry(xd, 'model')
            call pm%eval(pot)
            fd(:, n) = (pp%get_energy('adiabatic') - pm%get_energy('adiabatic')) / (2 * step)
            ! <a(x)|b(x + step)> - <a(x)|b(x - step)> = 2 step <a|d/dx b>.
            nac_fd(:, :, n) = (pp%link%s - pm%link%s) / (2 * step)
        end do

        ! Gradients.
        ret = 1
        err = 0.0_dp
        do a = 1, 4
            g = pot%get_gradient('adiabatic', a)
            err = max(err, maxval(abs(g - fd(a, :))))
        end do
        if (err > tol) return

        ! Nonadiabatic couplings within the states of the dynamics.
        ret = 2
        err = 0.0_dp
        do b = 1, 4
            do a = 1, 4
                if (a == b) cycle
                nac = pot%get_nadv('adiabatic', a, b)
                err = max(err, maxval(abs(nac - nac_fd(a, b, :)) / max(1.0_dp, abs(nac))))
            end do
        end do
        if (err > 1.0e-6_dp) return

        ! Hellmann-Feynman gradients miss the term from the change of the selected states.
        ret = 3
        call set_linear_model(pot_hf, h0, h1, x0, x1)
        call pot_hf%select_states([1, 1], [1, 3])
        pot_hf%exact_gradient = .false.
        call pot_hf%update_geometry(x, 'model')
        call pot_hf%eval()
        if (any(pot_hf%get_energy('adiabatic') /= pot%get_energy('adiabatic'))) return
        err_hf = 0.0_dp
        do a = 1, 4
            g = pot_hf%get_gradient('adiabatic', a)
            err_hf = max(err_hf, maxval(abs(g - fd(a, :))))
        end do
        if (err_hf < 100 * tol) return
        ret = 0
    end function subset_gradient_test


    !> Effective spin states with a subset of the states (S1 and T1 of the linear model). The
    !! effective model is built from all states of the source and takes the subset and the
    !! gradient option of its source. For couplings of rank 1, the adiabatic energies of the
    !! components of the subset are those of the projection model and the dark states at the
    !! energy of T1, twice.
    integer function effective_subset_test() result(ret)
        use global_defs, only : dp
        use pot_test_mod, only : linear_evaluator
        use evaluator_effective_spin_mod, only : effective_spin_evaluator
        type(linear_evaluator) :: src
        type(effective_spin_evaluator) :: eff
        real(dp) :: h0(3, 2, 2), h1(3, 2, 2, 2), x0(3, 3), x1(3, 3, 2)
        real(dp) :: geom(1, 2), u0(3), c(2, 2)
        real(dp), allocatable :: e_mch(:), e_src(:), e_eff(:)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: i, k

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
        call src%select_states([1, 1], [1, 3])
        src%exact_gradient = .false.
        call eff%init(src, 'projection')

        ret = 1
        if (eff%basis_size('diabatic') /= 4 .or. eff%basis_size('adiabatic') /= 2) return
        if (eff%basis_size('mch_sf') /= 2 .or. eff%exact_gradient) return
        if (any(shape(eff%proj) /= [8, 4]) .or. any(abs(eff%u - u0) > tol)) return

        ret = 2
        geom(1, :) = [0.3_dp, -0.2_dp]
        call src%update_geometry(geom, 'model')
        call src%eval()
        call eff%update_geometry(geom, 'model')
        call eff%eval()
        e_mch = src%get_energy('mch_sf')
        e_src = src%get_energy('adiabatic')
        e_eff = eff%get_energy('adiabatic')
        if (size(e_src) /= 4 .or. size(e_eff) /= 2) return
        if (any(abs(eff%get_energy('mch_sf') - e_mch) > tol)) return
        if (any(abs(sorted([e_eff, e_mch(2), e_mch(2)]) - e_src) > tol)) return
        if (any(shape(eff%get_transformation('diabatic', 'adiabatic')) /= [4, 2])) return
        ret = 0
    end function effective_subset_test


    !> Cases which stop with an error: in effective_init, a doublet, a source of effective states,
    !! no singlet-triplet coupling at the reference geometry and the projection model without
    !! signs; in model_init, an unknown parameter; in select_states, more states than the model
    !! has and a multiplicity which the model does not have.
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
            call src%set_model_states([1, 1], [1, 2], 'components')
            src%classical_coordinates = 'model'
        case('effective')
            call src%set_model_states([1, 1], [1, 3], 'effective_rss')
            src%classical_coordinates = 'model'
        case('no_coupling')
            h0 = 0.0_dp
            h1 = 0.0_dp
            h1(1, 1, 1, 1) = 1.0e-4_dp
            x0 = 0.0_dp
            x1 = 0.0_dp
            call set_linear_model(src, h0, h1, x0, x1)
        case('subset_too_many')
            call model%init([string('granucci-sst')])
            call model%select_states([3, 1], [1, 3])
            return
        case('subset_multiplicity')
            call model%init([string('granucci-sst')])
            call model%select_states([1, 1], [1, 2])
            return
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


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: external_json_test
    !> @brief The external_json evaluator with scripted responses of a program.
    !> @details
    !! Three states of two atoms. The program is 'cp' of a response written by the test.
    !! - Point A: energies, the requested gradient, wave functions, no link.
    !! - Point B (a copy of A): state 2 returned with the opposite sign. The sign is matched to
    !!   A, so link%s has a positive diagonal and the coupling of 1 and 2 changes sign. The
    !!   adiabatic, mch and mch_sf bases are the same states, and there is no diabatic basis.
    !! - complete on B: a missing gradient, and the coupling of 2 and 3 from a recalculation in
    !!   which state 3 has the opposite sign; the stored states keep their signs and the coupling
    !!   gets the sign of the stored states.
    !! - save_state/load_state of B into a new evaluator.
    !----------------------------------------------------------------------------------------------
    integer function external_json_test() result(ret)
        use global_defs, only : dp
        use evaluator_base_mod, only : new_request
        use external_json_mod, only : external_json_evaluator
        use state_overlap_mod, only : last_requested_id
        use json_module, only : json_core, json_value
        type(external_json_evaluator) :: a, b, c
        real(dp) :: geom(3, 2), g1(3, 2), g2(3, 2), v12(3, 2), v23(3, 2), r(3, 3), raw(3, 3)
        real(dp) :: e(3)
        type(json_core) :: json
        type(json_value), pointer :: root
        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: th

        geom = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp], [3, 2])
        g1 = reshape([0.1_dp, 0.2_dp, 0.3_dp, -0.1_dp, -0.2_dp, -0.3_dp], [3, 2])
        g2 = 2.0_dp * g1 + 0.05_dp
        v12 = reshape([1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, -2.0_dp], [3, 2])
        v23 = 3.0_dp * v12 - 0.1_dp
        e = [-1.0_dp, -0.8_dp, -0.7_dp]
        call setup_external(a, 'cp ext_scripted.json qm_out.json')

        ! Point A.
        ret = 1
        call write_scripted(last_requested_id() + 1, e, 1, g1, reshape([0], [2, 0]), v12,       &
        &                   unit3())
        call a%update_geometry(geom, 'cartesian')
        call a%set_request(new_request('adiabatic', gradient=[1]))
        call a%eval()
        if (a%link%valid()) return
        if (any(abs(a%get_energy('adiabatic') - e) > tol)) return
        if (any(abs(a%get_gradient('adiabatic', 1) - reshape(g1, [6])) > tol)) return
        if (a%have_grad(2) .or. any(a%have_nadv)) return
        if (any(abs(a%get_Hamiltonian('mch') - diag3(e)) > tol)) return
        if (any(a%get_energy('mch_sf') /= a%get_energy('adiabatic'))) return
        if (any(a%get_transformation('mch', 'adiabatic') /= unit3())) return
        if (any(a%transform('mch', 'adiabatic', diag3(e)) /= diag3(e))) return
        if (a%previous_point() /= 0 .or. a%id_diab /= 0 .or. a%diab%valid()) return

        ! Point B: a small rotation of the states, with state 2 of the opposite sign.
        ret = 2
        th = 0.05_dp
        r = unit3()
        r(1, 1) = cos(th)
        r(2, 1) = sin(th)
        r(1, 2) = -sin(th)
        r(2, 2) = cos(th)
        raw = r
        raw(:, 2) = -raw(:, 2)
        b = a
        call write_scripted(last_requested_id() + 1, e + 0.01_dp, 2, g2,                         &
        &                   reshape([1, 2], [2, 1]), v12, raw)
        call b%update_geometry(geom + 0.01_dp, 'cartesian')
        call b%set_request(new_request('adiabatic', gradient=[2], nadv=reshape([1, 2], [2, 1])))
        call b%eval(a)
        if (.not. b%link%matches(a%id, 'adiabatic', b%id, 'adiabatic')) return
        if (any(abs(b%link%s - r) > tol)) return
        if (any(b%get_transformation('adiabatic', 'adiabatic', from_point=a%id) /= b%link%t))    &
        &   return
        if (any(abs(b%get_transformation('adiabatic', 'mch', from_point=a%id) - r) > tol)) return
        if (any(abs(b%wf - r) > tol)) return
        if (any(abs(b%get_nadv('adiabatic', 1, 2) + reshape(v12, [6])) > tol)) return
        if (any(abs(b%get_nadv('adiabatic', 2, 1) - reshape(v12, [6])) > tol)) return
        if (b%have_grad(1) .or. .not. b%have_grad(2)) return

        ! complete: gradient of state 1, and the coupling of states 2 and 3 from a recalculation
        ! in which state 3 has the opposite sign.
        ret = 3
        raw = r
        raw(:, 3) = -raw(:, 3)
        call write_scripted(b%id, e + 0.01_dp, 1, g1, reshape([2, 3], [2, 1]), v23, raw)
        call b%complete(new_request('adiabatic', gradient=[1, 2], nadv=reshape([2, 3], [2, 1])))
        if (any(abs(b%wf - r) > tol)) return
        if (any(abs(b%get_gradient('adiabatic', 1) - reshape(g1, [6])) > tol)) return
        if (any(abs(b%get_gradient('adiabatic', 2) - reshape(g2, [6])) > tol)) return
        if (any(abs(b%get_nadv('adiabatic', 2, 3) + reshape(v23, [6])) > tol)) return
        if (any(abs(b%get_nadv('adiabatic', 1, 2) + reshape(v12, [6])) > tol)) return
        ! Nothing missing: no calculation (the scripted response is for another point).
        call write_scripted(b%id + 100, e, 1, g1, reshape([0], [2, 0]), v12, r)
        call b%complete(new_request('mch_sf', gradient=[2], nadv=reshape([3, 2], [2, 1])))

        ! Save and load.
        ret = 4
        call json%initialize()
        call json%create_object(root, '')
        call b%save_state(json, root)
        call setup_external(c, 'false')
        call c%load_state(json, root)
        call json%destroy(root)
        if (c%id /= b%id) return
        if (.not. c%link%matches(a%id, 'adiabatic', b%id, 'adiabatic')) return
        if (any(c%link%t /= b%link%t) .or. any(c%wf /= b%wf)) return
        if (any(c%adiab_vec /= b%adiab_vec)) return
        if (any(c%energy /= b%energy) .or. any(c%grad /= b%grad) .or. any(c%nadv /= b%nadv)) return
        if (any(c%have_grad .neqv. b%have_grad) .or. any(c%have_nadv .neqv. b%have_nadv)) return
        if (any(c%geometry /= b%geometry) .or. c%coordinates /= b%coordinates) return
        ret = 0
    end function external_json_test


    !> External evaluator of three singlets of two atoms, run with command.
    subroutine setup_external(ext, command)
        use external_json_mod, only : external_json_evaluator
        type(external_json_evaluator), intent(inout) :: ext
        character(len=*), intent(in) :: command

        call ext%init(command)
        call ext%set_model_states([3], [1], 'components')
        call ext%set_atoms([character(len=2) :: 'H', 'H'])
    end subroutine setup_external


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_scripted
    !> @brief Write a response for point id to ext_scripted.json, with the energies, the gradient
    !! of one state (none if grad_state is 0), the coupling vector of each pair of nacv_pair (all
    !! equal to nacv) and the wave functions wf.
    !----------------------------------------------------------------------------------------------
    subroutine write_scripted(id, energy, grad_state, grad, nacv_pair, nacv, wf)
        use global_defs, only : dp
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        integer, intent(in) :: id
        real(dp), intent(in) :: energy(:)
        integer, intent(in) :: grad_state
        real(dp), intent(in) :: grad(:, :) !< (3, natom)
        integer, intent(in) :: nacv_pair(:, :)
        real(dp), intent(in) :: nacv(:, :) !< (3, natom)
        real(dp), intent(in) :: wf(:, :)
        type(json_core) :: json
        type(json_value), pointer :: root, arr, obj
        integer :: k

        call json%initialize()
        call json%create_object(root, '')
        call json%add(root, 'protocol', 'zaghop')
        call json%add(root, 'version', 1)
        call json%add(root, 'point_id', id)
        call json%create_array(arr, 'states')
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'energies', energy)
        call json%add(arr, obj)
        call json%add(root, arr)
        call json%create_array(arr, 'gradients')
        if (grad_state > 0) then
            call json%create_object(obj, '')
            call json%add(obj, 'multiplicity', 1)
            call json%add(obj, 'state', grad_state)
            call json_add_array(json, obj, 'gradient', transpose(grad))
            call json%add(arr, obj)
        end if
        call json%add(root, arr)
        call json%create_array(arr, 'nacv')
        do k = 1, size(nacv_pair, 2)
            call json%create_object(obj, '')
            call json%add(obj, 'multiplicity', 1)
            call json%add(obj, 'states', nacv_pair(:, k))
            call json_add_array(json, obj, 'vector', transpose(nacv))
            call json%add(arr, obj)
        end do
        call json%add(root, arr)
        call json%create_array(arr, 'wave_function')
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'format', 'vectors')
        call json_add_array(json, obj, 'c', wf)
        call json%add(arr, obj)
        call json%add(root, arr)
        call json%print(root, 'ext_scripted.json')
        call json%destroy(root)
    end subroutine write_scripted


    pure function unit3() result(u)
        use global_defs, only : dp
        real(dp) :: u(3, 3)
        integer :: i

        u = 0.0_dp
        do i = 1, 3
            u(i, i) = 1.0_dp
        end do
    end function unit3


    pure function diag3(v) result(d)
        use global_defs, only : dp
        real(dp), intent(in) :: v(3)
        real(dp) :: d(3, 3)
        integer :: i

        d = 0.0_dp
        do i = 1, 3
            d(i, i) = v(i)
        end do
    end function diag3


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_error
    !> @brief Cases of the external_json evaluator which should stop with an error.
    !----------------------------------------------------------------------------------------------
    subroutine external_error(case)
        use global_defs, only : dp
        use evaluator_base_mod, only : new_request
        use external_json_mod, only : external_json_evaluator
        use state_overlap_mod, only : last_requested_id
        character(len=*), intent(in) :: case
        type(external_json_evaluator) :: ext
        real(dp) :: geom(3, 2), g(3, 2)
        real(dp), allocatable :: grad(:)

        geom = 0.0_dp
        geom(3, 2) = 1.4_dp
        g = 0.1_dp
        select case(case)
        case('external_no_response')
            call setup_external(ext, 'true')
        case('external_exit')
            call setup_external(ext, 'false')
        case('external_point_id')
            call setup_external(ext, 'cp ext_scripted.json qm_out.json')
            call write_scripted(last_requested_id() + 2, [0.0_dp, 0.1_dp, 0.2_dp], 1, g,         &
            &                   reshape([0], [2, 0]), g, unit3())
        case('external_missing_gradient')
            call setup_external(ext, 'cp ext_scripted.json qm_out.json')
            call write_scripted(last_requested_id() + 1, [0.0_dp, 0.1_dp, 0.2_dp], 2, g,         &
            &                   reshape([0], [2, 0]), g, unit3())
        case('external_not_computed')
            call setup_external(ext, 'cp ext_scripted.json qm_out.json')
            call write_scripted(last_requested_id() + 1, [0.0_dp, 0.1_dp, 0.2_dp], 1, g,         &
            &                   reshape([0], [2, 0]), g, unit3())
        case default
            write(*, *) 'Unknown case: ', case
            return
        end select
        call ext%update_geometry(geom, 'cartesian')
        call ext%set_request(new_request('adiabatic', gradient=[1]))
        call ext%eval()
        grad = ext%get_gradient('adiabatic', 2)
    end subroutine external_error


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

        call pot%set_model_states([2, 2], [1, 3], 'components')
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
