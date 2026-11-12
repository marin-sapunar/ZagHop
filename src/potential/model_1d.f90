!--------------------------------------------------------------------------------------------------
! MODULE: model_1d_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date September, 2024
!
! DESCRIPTION:
!> @brief One-dimensional analytic models for testing the code.
!> @details
!! The model is chosen with the first argument, followed by pairs of a parameter name and its
!! value (qlib model <name> [<parameter> <value> ...]):
!! - tully-i, tully-ii, tully-iii - Tully's models of two coupled singlets (J. Chem. Phys. 93,
!!   1061 (1990)), with parameters a, b, c, d (and e0 for tully-ii).
!! - granucci-st - One singlet and one triplet, whose spin-orbit couplings change sign at
!!   r_sigma (Granucci, Persico and Spighi, J. Chem. Phys. 137, 22A501 (2012), Sec. III A).
!!   Parameters r_sigma (default 10 bohr) and dr_sigma (width of the sign change, 2 bohr).
!! - granucci-sst - Two singlets and one triplet, with constant couplings of equal strength to
!!   the triplet, at the angle phi in the space of the triplet components (Sec. III B of the same
!!   article). Parameters phi (degrees, default 45) and de (E_S2 - E_S1, 0.035 Eh).

!! The triplet components of the Granucci models are those of the real basis of Eq. 12,
!!     T+ = (T(1) + T(-1)) / sqrt(2),  T- = i (T(1) - T(-1)) / sqrt(2),  Tz = i T(0),
!! in which the couplings z, ib and z* of a singlet with T(1), T(0) and T(-1) become the real
!! vector (sqrt(2) Re z, -sqrt(2) Im z, -b).
!--------------------------------------------------------------------------------------------------
module model_1d_mod
    use global_defs
    use string_mod
    use evaluator_diabatic_mod
    implicit none

    private
    public :: model_1d_evaluator


    ! Constants for Granucci models (Eqs. 17-19 of the article).
    real(dp), parameter :: g_a2 = 0.5_dp
    real(dp), parameter :: g_alpha1 = 0.35_dp
    real(dp), parameter :: g_alpha2 = 0.25_dp
    real(dp), parameter :: g_de = 0.04_dp
    real(dp), parameter :: g_rc = 10.0_dp
    !> a1 is set so that E_S crosses E_T at r_c (a1 = 0.03452 for r_c = 10).
    real(dp), parameter :: g_a1 = (g_a2 * exp(-g_alpha2 * g_rc) - g_de) * exp(g_alpha1 * g_rc)
    !> Couplings of the singlet with the triplet components in granucci-st before applying 
    !> the sign switch sigma(r). Eq. 11, Eq. 12 and Eq. 20 of the article.
    real(dp), parameter :: g_h_st(3) = [sqrt(2.0_dp) * 0.5e-3_dp, -sqrt(2.0_dp) * 0.5e-3_dp,     &
    &                                   -1.0e-3_dp]
    !> Strength of the couplings of each singlet with the triplet in granucci-sst.
    real(dp), parameter :: g_gamma_sst = sqrt(2.0_dp) * 1.0e-3_dp


    type, extends(diabatic_evaluator) :: model_1d_evaluator
        character(len=:), allocatable :: name
        real(dp), allocatable :: params(:)
        real(dp) :: x = 0.0_dp
    contains
        procedure :: init => model_init
        procedure :: set_geometry => model_set_geometry
        procedure :: eval_diab => model_eval_diab
        procedure :: default_reference => model_default_reference
        procedure, private :: read_params
    end type model_1d_evaluator

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: model_init
    !
    ! DESCRIPTION:
    !> @brief Initialize the model system with the given arguments.
    !> @details
    !! args(1) is the name of the model, followed by pairs of a parameter name and its value.
    !----------------------------------------------------------------------------------------------
    subroutine model_init(self, args)
        class(model_1d_evaluator) :: self
        type(string), intent(in) :: args(:)

        self%name = tolower(args(1)%s)
        self%classical_coordinates = 'model'
        self%ndof = 1

        select case(self%name)
        case('tully-i')
            call self%read_params(args(2:), [character(len=2) :: 'a', 'b', 'c', 'd'],             &
            &                     [0.01_dp, 1.6_dp, 0.005_dp, 1.0_dp])
        case('tully-ii')
            call self%read_params(args(2:), [character(len=2) :: 'a', 'b', 'c', 'd', 'e0'],       &
            &                     [0.1_dp, 0.28_dp, 0.015_dp, 0.06_dp, 0.05_dp])
        case('tully-iii')
            call self%read_params(args(2:), [character(len=2) :: 'a', 'b', 'c'],                  &
            &                     [6.0e-4_dp, 0.1_dp, 0.9_dp])
        case('granucci-st')
            call self%read_params(args(2:), [character(len=8) :: 'r_sigma', 'dr_sigma'],          &
            &                     [10.0_dp, 2.0_dp])
            if (self%params(2) <= 0.0_dp) call model_error('dr_sigma must be positive.')
        case('granucci-sst')
            call self%read_params(args(2:), [character(len=3) :: 'phi', 'de'], [45.0_dp, 0.035_dp])
        case default
            call model_error('Unrecognized model: '//self%name//'. Expected tully-i, '//          &
            &                'tully-ii, tully-iii, granucci-st or granucci-sst.')
        end select

        select case(self%name)
        case('granucci-st')
            ! One singlet and the three components of one triplet.
            call self%set_spin_manifolds([1, 1], [1, 3], 'components')
        case('granucci-sst')
            ! Two singlets and the three components of one triplet.
            call self%set_spin_manifolds([2, 1], [1, 3], 'components')
        case default
            ! Two singlet states, without spin-orbit coupling.
            call self%set_spin_manifolds([2], [1], 'components')
        end select
    end subroutine model_init


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_params
    !> @brief Set the parameters of the model from pairs of a name and a value.
    !> @details
    !! Parameters which are not given keep their default values.
    !----------------------------------------------------------------------------------------------
    subroutine read_params(self, args, names, defaults)
        class(model_1d_evaluator), intent(inout) :: self
        type(string), intent(in) :: args(:) !< Pairs of a parameter name and its value.
        character(len=*), intent(in) :: names(:) !< Names of the parameters.
        real(dp), intent(in) :: defaults(:) !< Default values of the parameters.
        integer :: i, j, ios
        character(len=:), allocatable :: expected

        self%params = defaults
        expected = trim(names(1))
        do j = 2, size(names)
            expected = expected//', '//trim(names(j))
        end do
        if (mod(size(args), 2) /= 0) then
            call model_error('The parameters of model '//self%name//' are given as pairs of a '//  &
            &                'name and a value.')
        end if
        do i = 1, size(args), 2
            j = findloc(names, tolower(args(i)%s), dim=1)
            if (j == 0) then
                call model_error('Unrecognized parameter of model '//self%name//': '//args(i)%s//  &
                &                '. Expected '//expected//'.')
            end if
            read(args(i+1)%s, *, iostat=ios) self%params(j)
            if (ios /= 0) then
                call model_error('Invalid value of parameter '//args(i)%s//': '//args(i+1)%s//'.')
            end if
        end do
    end subroutine read_params


    subroutine model_error(msg)
        character(len=*), intent(in) :: msg

        write(stderr, *) 'Error in model_1d_mod, model_init subroutine.'
        write(stderr, *) '  '//msg
        stop
    end subroutine model_error


    subroutine model_set_geometry(self, geometry, coordinates)
        class(model_1d_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        if (coordinates /= self%classical_coordinates) then
            write(stderr, *) 'Error in model_1d_mod, set_geometry subroutine.'
            write(stderr, *) '  Unsupported coordinates "'//coordinates//'". Expected "'//          &
            &                self%classical_coordinates//'".'
            stop
        end if

        self%x = geometry(1, 1)
    end subroutine model_set_geometry


    !> The default reference geometry is x = 0. For the Granucci models with the default
    !! parameters, all couplings there have their asymptotic values.
    subroutine model_default_reference(self, geom)
        class(model_1d_evaluator), intent(in) :: self
        real(dp), allocatable, intent(out) :: geom(:, :)

        allocate(geom(1, 1), source=0.0_dp)
    end subroutine model_default_reference


    subroutine model_eval_diab(self)
        class(model_1d_evaluator), intent(inout) :: self
        real(dp) :: v(3)
        real(dp) :: dv(3)

        select case(self%name)
        case('granucci-st')
            call granucci_st(self%params, self%x, self%sf_w, self%sf_dw(1, :, :), self%soc_w,      &
            &                self%soc_dw(1, :, :))
        case('granucci-sst')
            call granucci_sst(self%params, self%x, self%sf_w, self%sf_dw(1, :, :), self%soc_w,     &
            &                 self%soc_dw(1, :, :))
        case('tully-i')
            call tully_1(self%params, self%x, v, dv)
            call set_two_state(self, v, dv)
        case('tully-ii')
            call tully_2(self%params, self%x, v, dv)
            call set_two_state(self, v, dv)
        case('tully-iii')
            call tully_3(self%params, self%x, v, dv)
            call set_two_state(self, v, dv)
        end select
    end subroutine model_eval_diab


    !> Set the diabatic Hamiltonian and its gradient (1 mode) from V11, V22, V12 and derivatives.
    subroutine set_two_state(self, v, dv)
        class(model_1d_evaluator), intent(inout) :: self
        real(dp), intent(in) :: v(3) !< Vij matrix elements V11, V22, V12.
        real(dp), intent(in) :: dv(3) !< Derivative of Vij matrix elements.

        self%sf_w(1, 1) = v(1)
        self%sf_w(2, 2) = v(2)
        self%sf_w(1, 2) = v(3)
        self%sf_w(2, 1) = v(3)

        self%sf_dw(1, 1, 1) = dv(1)
        self%sf_dw(1, 2, 2) = dv(2)
        self%sf_dw(1, 1, 2) = dv(3)
        self%sf_dw(1, 2, 1) = dv(3)
    end subroutine set_two_state


    subroutine tully_1(p, x, v, dv)
        real(dp), intent(in) :: p(4) !< Model parameters A, B, C & D.
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: v(3) !< Vij matrix elements V11, V22, V12.
        real(dp), intent(out) :: dv(3) !< Derivative of Vij matrix elements.

        if (x >= 0.0_dp) then
            v(1) = p(1) * (1.0_dp - exp(-p(2) * x))
            dv(1) = p(2) * (p(1) - v(1))
        else
            v(1) = - p(1) * (1.0_dp - exp(p(2) * x))
            dv(1) = p(2) * (p(1) + v(1))
        end if
        v(2) = -v(1)
        v(3) = p(3) * exp(-p(4) * x**2.0_dp)
        dv(2) = - dv(1)
        dv(3) = - 2.0_dp * p(4) * x * v(3)
    end subroutine tully_1


    subroutine tully_2(p, x, v, dv)
        real(dp), intent(in) :: p(5) !< Model parameters A, B, C, D & E0.
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: v(3) !< Vij matrix elements V11, V22, V12.
        real(dp), intent(out) :: dv(3) !< Derivative of Vij matrix elements.

        v(1) = 0.0_dp
        v(2) = - p(1) * exp(- p(2) * x**2) + p(5)
        v(3) = p(3) * exp(- p(4) * x**2)
        dv(1) = 0.0_dp
        dv(2) = - 2.0_dp * p(2) * x * (v(2) - p(5))
        dv(3) = - 2.0_dp * p(4) * x * v(3)
    end subroutine tully_2


    subroutine tully_3(p, x, v, dv)
        real(dp), intent(in) :: p(3) !< Model parameters A, B & C.
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: v(3) !< Vij matrix elements V11, V22, V12.
        real(dp), intent(out) :: dv(3) !< Derivative of Vij matrix elements.

        v(1) = p(1)
        v(2) = -p(1)
        dv(1) = 0.0_dp
        dv(2) = 0.0_dp
        if (x >= 0.0_dp) then
            v(3) = p(2) * (2.0_dp - exp(- p(3) * x))
            dv(3) = p(3) * (2.0_dp * p(2) - v(3))
        else
            v(3) = p(2) * exp(p(3) * x)
            dv(3) = p(3) * v(3)
        end if
    end subroutine tully_3


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: granucci_energies
    !> @brief Singlet and triplet energies of the Granucci models and their derivatives.
    !----------------------------------------------------------------------------------------------
    pure subroutine granucci_energies(x, es, des, et, det)
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: es, des !< Singlet energy and its derivative.
        real(dp), intent(out) :: et, det !< Triplet energy and its derivative.

        es = g_a1 * exp(-g_alpha1 * x) + g_de
        des = -g_alpha1 * (es - g_de)
        et = g_a2 * exp(-g_alpha2 * x)
        det = -g_alpha2 * et
    end subroutine granucci_energies


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sign_switch
    !> @brief Function sigma(r) of Eq. 21 of the article, which changes smoothly from 1 to -1
    !! between r0 - dr/2 and r0 + dr/2, and its derivative.
    !----------------------------------------------------------------------------------------------
    pure subroutine sign_switch(x, r0, dr, s, ds)
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(in) :: r0 !< Zero of sigma.
        real(dp), intent(in) :: dr !< Width of the switch.
        real(dp), intent(out) :: s, ds
        real(dp) :: y

        y = (x - r0) / dr
        if (y <= -0.5_dp) then
            s = 1.0_dp
            ds = 0.0_dp
        else if (y >= 0.5_dp) then
            s = -1.0_dp
            ds = 0.0_dp
        else
            s = 4.0_dp * y**3 - 3.0_dp * y
            ds = (12.0_dp * y**2 - 3.0_dp) / dr
        end if
    end subroutine sign_switch


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: granucci_st
    !> @brief One singlet (state 1) and one triplet (states 2-4: T+, T-, Tz), with couplings
    !! sigma(r) h, where h is a fixed vector (Sec. III A of the article).
    !----------------------------------------------------------------------------------------------
    pure subroutine granucci_st(p, x, sf_w, sf_dw, soc_w, soc_dw)
        real(dp), intent(in) :: p(2) !< Parameters r_sigma and dr_sigma.
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: sf_w(2, 2), sf_dw(2, 2) !< Spin-free Hamiltonian, derivative.
        real(dp), intent(out) :: soc_w(4, 4), soc_dw(4, 4) !< Spin-orbit couplings, derivative.
        real(dp) :: s, ds

        sf_w = 0.0_dp
        sf_dw = 0.0_dp
        call granucci_energies(x, sf_w(1, 1), sf_dw(1, 1), sf_w(2, 2), sf_dw(2, 2))
        call sign_switch(x, p(1), p(2), s, ds)
        soc_w = 0.0_dp
        soc_dw = 0.0_dp
        soc_w(1, 2:4) = s * g_h_st
        soc_w(2:4, 1) = s * g_h_st
        soc_dw(1, 2:4) = ds * g_h_st
        soc_dw(2:4, 1) = ds * g_h_st
    end subroutine granucci_st


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: granucci_sst
    !> @brief Two singlets (states 1 and 2) and one triplet (states 3-5: T+, T-, Tz), with
    !! constant couplings (Sec. III B of the article).
    !> @details
    !! E_S2 = E_S, and E_S1 = E_S - de, parallel to it. The couplings of S1 and S2 with the
    !! triplet components are gamma (1, 0, 0) and gamma (cos phi, sin phi, 0), so phi is the angle
    !! phi_12 of Eq. 29 of the article. The singlets are not coupled to each other.
    !----------------------------------------------------------------------------------------------
    pure subroutine granucci_sst(p, x, sf_w, sf_dw, soc_w, soc_dw)
        real(dp), intent(in) :: p(2) !< Parameters phi (degrees) and de.
        real(dp), intent(in) :: x !< Position.
        real(dp), intent(out) :: sf_w(3, 3), sf_dw(3, 3) !< Spin-free Hamiltonian, derivative.
        real(dp), intent(out) :: soc_w(5, 5), soc_dw(5, 5) !< Spin-orbit couplings, derivative.
        real(dp) :: phi

        sf_w = 0.0_dp
        sf_dw = 0.0_dp
        call granucci_energies(x, sf_w(2, 2), sf_dw(2, 2), sf_w(3, 3), sf_dw(3, 3))
        sf_w(1, 1) = sf_w(2, 2) - p(2)
        sf_dw(1, 1) = sf_dw(2, 2)
        phi = p(1) * pi / 180.0_dp
        soc_w = 0.0_dp
        soc_dw = 0.0_dp
        soc_w(1, 3) = g_gamma_sst
        soc_w(2, 3:4) = g_gamma_sst * [cos(phi), sin(phi)]
        soc_w(3:5, 1:2) = transpose(soc_w(1:2, 3:5))
    end subroutine granucci_sst


end module model_1d_mod
