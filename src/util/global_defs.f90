!--------------------------------------------------------------------------------------------------
! MODULE: Global_defs
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date January, 2017
!
!> @brief Important parameters and commonly used simple types.
!--------------------------------------------------------------------------------------------------
module global_defs
    use, intrinsic :: iso_fortran_env, only : &
        stderr => error_unit, &
        stdout => output_unit, &
        stdin => input_unit, &
        sp => real32, & 
        dp => real64, & 
        qp => real128, &
        i4 => int32, &
        i8 => int64


    integer, parameter :: j15 = selected_int_kind(15)
    integer :: print_level !< Global print level option for programs.

    ! Number constants.
    real(dp), parameter :: num0 = 0.0_dp
    real(dp), parameter :: num1 = 1.0_dp
    real(dp), parameter :: num2 = 2.0_dp
    real(dp), parameter :: num3 = 3.0_dp
    real(dp), parameter :: num4 = 4.0_dp
    real(dp), parameter :: f1o2 = num1 / num2
    real(dp), parameter :: f1o4 = num1 / num4
    real(dp), parameter :: f2o3 = num2 / num3
    real(dp), parameter :: f3o4 = num3 / num4
    real(dp), parameter :: sqr2 = sqrt(num2)
    real(dp), parameter :: rsq2 = num1 / sqr2
    real(dp), parameter :: sqr3 = sqrt(num3)
    real(dp), parameter :: rsq3 = num1 / sqr3
    real(dp), parameter :: pi = acos(-num1)
    real(dp), parameter :: pio2 = pi / num2


 !  !----------------------------------------------------------------------------------------------
 !  ! TYPE: IVec
 !  !> @brief Array of integer vectors of varying dimensions.
 !  !----------------------------------------------------------------------------------------------
 !  type ivec 
 !     integer, allocatable :: c(:) !< Coefficents of the vectors.
 !  end type ivec


    !----------------------------------------------------------------------------------------------
    ! TYPE: IMat
    !> @brief Array of integer matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type imat
        integer, allocatable :: c(:, :) !< Coefficients of the matrices.
    end type imat


    !----------------------------------------------------------------------------------------------
    ! TYPE: RVec
    !> @brief Array of real vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type rvec 
       real(dp), allocatable :: c(:) !< Coefficents of the vectors.
    end type rvec


    !----------------------------------------------------------------------------------------------
    ! TYPE: RMat
    !> @brief Array of real matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type rmat
        real(dp), allocatable :: c(:, :) !< Coefficients of the matrices.
    end type rmat


    !----------------------------------------------------------------------------------------------
    ! TYPE: LVec
    !> @brief Array of logical vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type lvec
       logical, allocatable :: l(:)
    end type lvec


    !----------------------------------------------------------------------------------------------
    ! TYPE: LMat
    !> @brief Array of logical matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type lmat
       logical, allocatable :: l(:, :)
    end type lmat


contains


    subroutine errstop(errsub, errmsg, errnum)
        character(len=*), intent(in) :: errsub !< Subroutine where the error occured.
        character(len=*), intent(in), optional :: errmsg !< Error message.
        integer, intent(in), optional :: errnum !< Error number.

        write(stderr, '(1x,a)') 'Error in '//errsub//'.'
        if (present(errmsg)) write(stderr, '(3x, a, a)')  'Error: ', errmsg
        if (present(errnum)) write(stderr, '(3x, a, i0)')  'Status: ', errnum
        flush(stderr)
        call abort()
    end subroutine errstop


    elemental function truefalse_str(bool) result(tf)
        logical, intent(in) :: bool
        character(len=5) :: tf

        if (bool) then
            tf = "True "
        else
            tf = "False"
        end if
    end function truefalse_str


end module global_defs
