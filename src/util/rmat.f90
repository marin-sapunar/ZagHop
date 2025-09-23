module rmat_mod
    use global_defs
    implicit none

    private
    public :: rmat

    type rmat
        real(dp), allocatable :: c(:, :) !< Matrix values.
    contains
        procedure :: size => size_rmat
    end type rmat

    public :: assignment(=)
    interface assignment(=)
        module procedure assign_rmat_rmat
        module procedure assign_rmat_array
        module procedure assign_array_rmat
    end interface assignment(=)

contains

    subroutine assign_rmat_rmat(mat1, mat2)
        class(rmat), intent(inout) :: mat1
        class(rmat), intent(in) :: mat2

        if (allocated(mat1%c)) deallocate(mat1%c)
        allocate(mat1%c(size(mat2%c, 1), size(mat2%c, 2)))
        mat1%c(:, :) = mat2%c(:, :)
    end subroutine assign_rmat_rmat

    subroutine assign_rmat_array(rm, arr)
        class(rmat), intent(out) :: rm
        real(dp), intent(in) :: arr(:, :)

        if (allocated(rm%c)) deallocate(rm%c)
        allocate(rm%c(size(arr, 1), size(arr, 2)))
        rm%c(:, :) = arr(:, :)
    end subroutine assign_rmat_array

    subroutine assign_array_rmat(arr, mat)
        real(dp), intent(out) :: arr(:, :)
        class(rmat), intent(in) :: mat

        arr = mat%c
    end subroutine assign_array_rmat

    function size_rmat(self, dim) result(sz)
        class(rmat), intent(in) :: self
        integer, intent(in), optional :: dim
        integer :: sz

        if (present(dim)) then
            if (dim < 1 .or. dim > 2) then
                write(*, *) "Error in rmat size: dimension must be 1 or 2."
                stop
            end if
            if (.not. allocated(self%c)) then
                sz = 0
            else
                sz = size(self%c, dim)
            end if
        else
            if (.not. allocated(self%c)) then
                sz = 0
            else
                sz = size(self%c)
            end if
        end if
    end function size_rmat

end module rmat_mod