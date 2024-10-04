!--------------------------------------------------------------------------------------------------
! MODULE: model_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date September, 2024
!
! DESCRIPTION: 
!> @brief Simple analytic models for testing the code.
!--------------------------------------------------------------------------------------------------
module model_mod
    use global_defs
    use string_mod
    implicit none

    private
    public :: qmodel
    public :: model_system


    type model_system
        character(len=:), allocatable :: name
        real(dp), allocatable :: params(:)
    contains
        procedure :: init => model_init
        procedure :: eval => model_eval_sym2d
    end type model_system

    type(model_system) :: qmodel

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: model_init
    !
    ! DESCRIPTION:
    !> @brief
    !> @details
    !----------------------------------------------------------------------------------------------
    subroutine model_init(self, args)
        class(model_system) :: self
        type(string), intent(in) :: args(:)
        integer :: narg
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
    end subroutine model_init


    subroutine model_eval_sym2d(self, x, cstate, en, grad, nadv)
        class(model_system) :: self
        real(dp), intent(in) :: x(1, 1)
        integer, intent(in) :: cstate
        real(dp), intent(out) :: en(2)
        real(dp), intent(out) :: grad(1, 1)
        real(dp), allocatable, intent(inout) :: nadv(:, :, :)
        real(dp) :: v(3) !< Vij matrix elements V11, V22, V12.
        real(dp) :: dv(3) !< Derivative of Vij matrix elements.
        real(dp) :: diag_dif
        real(dp) :: diag_dif_g
        real(dp) :: sq
        real(dp) :: sum_gr

        select case(self%name)
        case('tully-i')
            call tully_1(self%params, x(1, 1), v, dv)
        case('tully-ii')
            call tully_2(self%params, x(1, 1), v, dv)
        case('tully-iii')
            call tully_3(self%params, x(1, 1), v, dv)
        end select
 
        diag_dif = (v(2) - v(1)) * 0.5_dp
        diag_dif_g = (dv(2) - dv(1)) * 0.5_dp
        sq = sqrt(diag_dif**2 + v(3)**2)
        sum_gr = diag_dif * diag_dif_g + v(3) * dv(3)

        en(1) = (v(1) + v(2)) * 0.5_dp - sq
        en(2) = (v(1) + v(2)) * 0.5_dp + sq

        if (cstate == 1) then
            grad(1, 1) = (dv(1) + dv(2)) * 0.5_dp - sum_gr / sq
        else
            grad(1, 1) = (dv(1) + dv(2)) * 0.5_dp + sum_gr / sq
        end if

        if (allocated(nadv)) then
            nadv = 0.0_dp
            nadv(1, 1, 2) = 0.5_dp / (1 + v(3)**2 / diag_dif**2) 
            nadv(1, 1, 2) = nadv(1, 1, 2) * (dv(3) / diag_dif - v(3) * diag_dif_g / diag_dif**2)
            nadv(1, 2, 1) = - nadv(1, 1, 2)
        end if
    end subroutine model_eval_sym2d


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
    ! SUBROUTINE: model_tully_1
    !
    ! DESCRIPTION:
    !> @brief
    !> @details
    !! 
    !----------------------------------------------------------------------------------------------
  ! subroutine model_tully_1(x, cstate, en, grad, nadv)
  !     real(dp), intent(in) :: x !< Position.
  !     integer, intent(in) :: cstate !< Currently populated state.
  !     real(dp), intent(out) :: en(2) !< Adiabatic energies.
  !   ! real(dp), intent(out) :: grad(1, 2) !< Adiabatic gradients.
  !     real(dp), intent(out) :: grad(1, 1) !< Adiabatic gradients.
  !     real(dp), intent(out) :: nadv(1, 2, 2) !< Nonadiabatic coupling vectors.
  !     real(dp), parameter :: a = 0.01_dp
  !     real(dp), parameter :: b = 1.6_dp
  !     real(dp), parameter :: c = 0.005_dp
  !     real(dp), parameter :: d = 1.0_dp
  !     real(dp) :: v11, v12, dv11, dv12
  !     real(dp) :: rnum(5)

  !     call random_number(rnum)
  !     rnum = (rnum - 0.5_dp) * 1.0e-7_dp

  !     if (x >= 0.0_dp) then
  !         v11 = a * (1.0_dp - exp(-b * x))
  !     else
  !         v11 = - a * (1.0_dp - exp(b * x))
  !     end if
  !     v12 = c * exp(-d * x**2.0_dp)
  !     dv11 = a * b * exp(-b * abs(x))
  !     dv12 = - 2.0_dp * d * x * v12

  !     en(2) = sqrt(v11**2 + v12**2)
  !     en(1) = - en(2)
  !    !grad(1, 2) = (v11 * dv11 + v12 * dv12) / en(2)
  !    !grad(1, 1) = -grad(1, 2)
  !     grad(1, 1) = (v11 * dv11 + v12 * dv12) / en(2)
  !     if (cstate == 1) grad(1, 1) = -grad(1, 1)
  !     en = en + rnum(1:2)
  !     grad = grad + rnum(3)
  !     nadv = 0.0_dp
  !     nadv(1, 1, 2) = (v12 * dv11 - v11 * dv12) * 0.5_dp / en(1)**2
  !     nadv(1, 2, 1) = - nadv(1, 1, 2)
  !     nadv(1, 1, 2) = nadv(1, 1, 2) + rnum(4)
  !     nadv(1, 2, 1) = nadv(1, 2, 1) + rnum(4)
  ! end subroutine model_tully_1


end module model_mod
