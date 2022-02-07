!-------------------------------------------------------------------------------------------------
! MODULE: ivec_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Array of integer vectors of varying dimensions.
!-------------------------------------------------------------------------------------------------
module ivec_mod
    implicit none

    private
    public :: ivec

    !----------------------------------------------------------------------------------------------
    ! TYPE: IVec
    !> @brief Array of integer vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type ivec 
       integer, allocatable :: c(:) !< Coefficents of the vectors.
    end type ivec


    public :: assignment(=)
    interface assignment(=)
        module procedure assign_iv_i
        module procedure assign_i_iv
    end interface assignment(=)


    public :: maxval
    interface maxval
        module procedure maxval_ivec
    end interface maxval


    public :: minval
    interface minval
        module procedure minval_ivec
    end interface minval


contains


    subroutine assign_iv_i(iv, i)
        type(ivec), intent(out) :: iv
        integer, intent(in) :: i(:)
        iv%c = i
    end subroutine assign_iv_i


    subroutine assign_i_iv(i, iv)
        integer, intent(out) :: i(:)
        type(ivec), intent(in) :: iv
        i = iv%c
    end subroutine assign_i_iv


    function maxval_ivec(iv) result(val)
        type(ivec), intent(in) :: iv(:)
        integer :: val
        integer :: i, mxi

        val = iv(1)%c(1)
        do i = 1, size(iv)
            mxi = maxval(iv(i)%c)
            if (mxi > val) val = mxi
        end do
    end function maxval_ivec


    function minval_ivec(iv) result(val)
        type(ivec), intent(in) :: iv(:)
        integer :: val
        integer :: i, mxi

        val = iv(1)%c(1)
        do i = 1, size(iv)
            mxi = minval(iv(i)%c)
            if (mxi < val) val = mxi
        end do
    end function minval_ivec


end module ivec_mod
