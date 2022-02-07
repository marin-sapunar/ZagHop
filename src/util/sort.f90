!-------------------------------------------------------------------------------------------------
! MODULE: sort_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Basic (slow) sorting routines.
!-------------------------------------------------------------------------------------------------
module sort_mod
    use global_defs
    implicit none


    private
    public :: sort

    interface sort
        module procedure sort_selection_int
        module procedure sort_selection_rel
    end interface sort

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sort_selection_int
    !> @brief In place selection sort of a list of integers
    !----------------------------------------------------------------------------------------------
    subroutine sort_selection_int(a)
        integer, intent(inout) :: a(:) !< List to be sorted.
        integer :: pos
        integer :: i
        integer :: temp

        do i = 1, size(a)-1
            pos = minloc(a(i:), 1) + i - 1
            if (a(i) > a(pos)) then
                temp = a(i)
                a(i) = a(pos)
                a(pos) = temp
            end if
        end do
    end subroutine sort_selection_int
   

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sort_selection_rel
    !> @brief In place selection sort of a list of reals.
    !----------------------------------------------------------------------------------------------
    subroutine sort_selection_rel(a)
        real(dp), intent(inout) :: a(:) !< List to be sorted.
        integer :: pos
        integer :: i
        real(dp) :: temp

        do i = 1, size(a)-1
            pos = minloc(a(i:), 1) + i - 1
            if (a(i) > a(pos)) then
                temp = a(i)
                a(i) = a(pos)
                a(pos) = temp
            end if
        end do
    end subroutine sort_selection_rel


end module sort_mod
