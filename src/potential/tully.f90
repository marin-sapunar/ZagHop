 !--------------------------------------------------------------------------------------------------
! MODULE: tully_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date September, 2024
!
! DESCRIPTION: 
!> @brief Simple analytic models for testing the code.
!--------------------------------------------------------------------------------------------------
module tully_mod
    use global_defs
    use string_mod
    use evaluator_diabatic_mod
    implicit none

    private
    public :: tully_model


    type, extends(diabatic_evaluator) :: tully_model
        character(len=:), allocatable :: name
        real(dp), allocatable :: params(:)
        real(dp) :: x = 0.0_dp
        real(dp), allocatable :: saved_nadv(:, :, :)
    contains
        procedure :: init => model_init
        procedure :: update_geometry => model_update_geometry
        procedure :: eval => model_eval
        procedure :: get_nadv => model_get_nadv
    end type tully_model

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: model_init
    !
    ! DESCRIPTION:
    !> @brief Initialize the model system with the given arguments.
    !> @details
    !----------------------------------------------------------------------------------------------
    subroutine model_init(self, args)
        class(tully_model) :: self
        type(string), intent(in) :: args(:)
        integer :: i

        self%name = tolower(args(1)%s)

        select case(self%name)
        case('tully-i')
            allocate(self%params(4))
            self%params(1) = 0.01_dp ! A
            self%params(2) = 1.6_dp ! B
            self%params(3) = 0.005_dp ! C
            self%params(4) = 1.0_dp ! D
            do i = 2, size(args), 2
                select case(args(i)%s)
                case('a','A')
                    read(args(i+1)%s, *) self%params(1)
                case('b','B')
                    read(args(i+1)%s, *) self%params(2)
                case('c','C')
                    read(args(i+1)%s, *) self%params(3)
                case('d','D')
                    read(args(i+1)%s, *) self%params(4)
                end select
            end do
        case('tully-ii')
            allocate(self%params(5))
            self%params(1) = 0.1_dp ! A
            self%params(2) = 0.28_dp ! B
            self%params(3) = 0.015_dp ! C
            self%params(4) = 0.06_dp ! D
            self%params(5) = 0.05_dp ! E0
            do i = 2, size(args), 2
                select case(args(i)%s)
                case('a','A')
                    read(args(i+1)%s, *) self%params(1)
                case('b','B')
                    read(args(i+1)%s, *) self%params(2)
                case('c','C')
                    read(args(i+1)%s, *) self%params(3)
                case('d','D')
                    read(args(i+1)%s, *) self%params(4)
                case('e0','E0')
                    read(args(i+1)%s, *) self%params(5)
                end select
            end do
        case('tully-iii')
            allocate(self%params(3))
            self%params(1) = 6.0e-4_dp ! A
            self%params(2) = 0.1_dp ! B
            self%params(3) = 0.9_dp ! C
            do i = 2, size(args), 2
                select case(args(i)%s)
                case('a','A')
                    read(args(i+1)%s, *) self%params(1)
                case('b','B')
                    read(args(i+1)%s, *) self%params(2)
                case('c','C')
                    read(args(i+1)%s, *) self%params(3)
                end select
            end do
        end select

        ! Initialize diabatic_evaluator arrays for 2 states, 1 group.
        self%n_state = 2
        self%n_group = 1
        allocate(self%group_nstate(1), source=2)
        allocate(self%group_i0(1), source=1)
        allocate(self%diab_w(2, 2))
        allocate(self%diab_dw(1, 2, 2))
        allocate(self%adiab_w(2, 2))
        allocate(self%adiab_dw(1, 2, 2))
        allocate(self%group_adiab_w(2, 2))
        allocate(self%group_adiab_dw(1, 2, 2))
        allocate(self%group_adiab_trans(2, 2))
        allocate(self%saved_nadv(1, 2, 2))
    end subroutine model_init


    subroutine model_update_geometry(self, geometry)
        class(tully_model), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)

        self%x = geometry(1, 1)
        if (allocated(self%adiab_trans)) self%ldiab_trans = self%adiab_trans
    end subroutine model_update_geometry


    subroutine model_eval(self)
        class(tully_model), intent(inout) :: self
        real(dp) :: v(3)
        real(dp) :: dv(3)
        real(dp) :: edif
        integer :: imode
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        ! Compute diabatic matrix elements.
        select case(self%name)
        case('tully-i')
            call tully_1(self%params, self%x, v, dv)
        case('tully-ii')
            call tully_2(self%params, self%x, v, dv)
        case('tully-iii')
            call tully_3(self%params, self%x, v, dv)
        end select

        ! Set diabatic Hamiltonian.
        self%diab_w(1, 1) = v(1)
        self%diab_w(2, 2) = v(2)
        self%diab_w(1, 2) = v(3)
        self%diab_w(2, 1) = v(3)

        ! Set diabatic gradient (1 mode).
        self%diab_dw(1, 1, 1) = dv(1)
        self%diab_dw(1, 2, 2) = dv(2)
        self%diab_dw(1, 1, 2) = dv(3)
        self%diab_dw(1, 2, 1) = dv(3)

        ! Diagonalize to get adiabatic transformation.
        call self%eval_eigvec('adiabatic')
        call self%eval_eigvec('group_adiabatic')

        ! Phase matching with previous step.
        if (allocated(self%ldiab_trans)) then
            call match_phase(self%ldiab_trans, self%adiab_trans)
        end if

        ! Transform to adiabatic basis.
        self%adiab_w = self%diab_w
        self%adiab_dw = self%diab_dw
        call self%change_basis('adiabatic', self%adiab_w)
        do imode = 1, 1
            call self%change_basis('adiabatic', self%adiab_dw(imode, :, :))
        end do

        self%group_adiab_w = self%diab_w
        self%group_adiab_dw = self%diab_dw
        call self%change_basis('group_adiabatic', self%group_adiab_w)
        do imode = 1, 1
            call self%change_basis('group_adiabatic', self%group_adiab_dw(imode, :, :))
        end do

        ! Compute and store nonadiabatic coupling vectors.
        self%saved_nadv = 0.0_dp
        edif = self%adiab_w(2, 2) - self%adiab_w(1, 1)
        if (abs(edif) < tiny_hf) then
            edif = sign(tiny_hf, edif)
        end if
        self%saved_nadv(1, 1, 2) = self%adiab_dw(1, 1, 2) / edif
        self%saved_nadv(1, 2, 1) = -self%saved_nadv(1, 1, 2)
    end subroutine model_eval


    function model_get_nadv(self, basis, istate1, istate2) result(nadv)
        class(tully_model), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)

        allocate(nadv(1))
        nadv(1) = self%saved_nadv(1, istate1, istate2)
    end function model_get_nadv


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


end module tully_mod
