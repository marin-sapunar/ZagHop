!--------------------------------------------------------------------------------------------------
! MODULE: json_util_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2026
!
! DESCRIPTION:
!> @brief Helper procedures for reading and writing JSON values with json-fortran.
!> @details
!! json-fortran (9.3.1) reads and writes scalars and vectors, but not arrays of higher rank or
!! complex numbers. Here, a rank 2 or 3 array is a nested JSON array with the Fortran shape:
!! element a(i, j, k) is a[i][j][k], so numpy.array() of the JSON value has the shape of the
!! Fortran array. A complex vector is an object with the real and imaginary parts, {"re": [...],
!! "im": [...]}.
!!
!! The json_get_required procedures stop the program with the path of the value if it is missing
!! or has the wrong type or shape. They are meant for files written by the program itself, such
!! as backup files, where a missing value is an error.
!--------------------------------------------------------------------------------------------------
module json_util_mod
    use global_defs
    use json_module, only : json_core, json_value, json_array
    implicit none

    private
    public :: json_add_array
    public :: json_get_required
    public :: json_get_child
    public :: json_check_value


    interface json_add_array
        module procedure add_real_2d
        module procedure add_real_3d
        module procedure add_complex_1d
    end interface json_add_array


    interface json_get_required
        module procedure get_integer
        module procedure get_real
        module procedure get_logical
        module procedure get_string
        module procedure get_integer_1d
        module procedure get_real_1d
        module procedure get_real_2d
        module procedure get_real_3d
        module procedure get_complex_1d
    end interface json_get_required


    interface json_check_value
        module procedure check_integer
        module procedure check_string
        module procedure check_logical
    end interface json_check_value


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: add_real_2d
    !> @brief Add a rank 2 real array to obj as a nested array, a[i][j] = a(i, j).
    !----------------------------------------------------------------------------------------------
    subroutine add_real_2d(json, obj, name, a)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the array is added.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a(:, :)
        type(json_value), pointer :: arr
        integer :: i

        call json%create_array(arr, name)
        do i = 1, size(a, 1)
            call json%add(arr, '', a(i, :))
        end do
        call json%add(obj, arr)
    end subroutine add_real_2d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: add_real_3d
    !> @brief Add a rank 3 real array to obj as a nested array, a[i][j][k] = a(i, j, k).
    !----------------------------------------------------------------------------------------------
    subroutine add_real_3d(json, obj, name, a)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the array is added.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a(:, :, :)
        type(json_value), pointer :: arr
        integer :: i

        call json%create_array(arr, name)
        do i = 1, size(a, 1)
            call add_real_2d(json, arr, '', a(i, :, :))
        end do
        call json%add(obj, arr)
    end subroutine add_real_3d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: add_complex_1d
    !> @brief Add a complex vector to obj as an object {"re": [...], "im": [...]}.
    !----------------------------------------------------------------------------------------------
    subroutine add_complex_1d(json, obj, name, a)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the vector is added.
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: a(:)
        type(json_value), pointer :: c

        call json%create_object(c, name)
        call json%add(c, 're', real(a, dp))
        call json%add(c, 'im', aimag(a))
        call json%add(obj, c)
    end subroutine add_complex_1d


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: json_get_child
    !> @brief Pointer to the value at path within obj. Stops if it is missing.
    !----------------------------------------------------------------------------------------------
    function json_get_child(json, obj, path) result(child)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        type(json_value), pointer :: child
        logical :: found

        nullify(child)
        call json%get(obj, path, child, found)
        if (.not. found .or. .not. associated(child)) call get_error(json, path, 'not found')
    end function json_get_child


    subroutine get_integer(json, obj, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        integer, intent(out) :: val
        logical :: found

        call json%get(obj, path, val, found)
        if (.not. found) call get_error(json, path, 'not found or not an integer')
    end subroutine get_integer


    subroutine get_real(json, obj, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: val
        logical :: found

        call json%get(obj, path, val, found)
        if (.not. found) call get_error(json, path, 'not found or not a number')
    end subroutine get_real


    subroutine get_logical(json, obj, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        logical, intent(out) :: val
        logical :: found

        call json%get(obj, path, val, found)
        if (.not. found) call get_error(json, path, 'not found or not a logical')
    end subroutine get_logical


    subroutine get_string(json, obj, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        character(len=:), allocatable, intent(out) :: val
        logical :: found

        call json%get(obj, path, val, found)
        if (.not. found) call get_error(json, path, 'not found or not a string')
    end subroutine get_string


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_integer_1d
    !> @brief Get an integer vector. With n, its size must be n.
    !----------------------------------------------------------------------------------------------
    subroutine get_integer_1d(json, obj, path, val, n)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        integer, allocatable, intent(out) :: val(:)
        integer, intent(in), optional :: n !< Expected size.
        type(json_value), pointer :: p

        p => json_get_child(json, obj, path)
        call check_array(json, p, path)
        call json%get(p, val)
        if (json%failed() .or. .not. allocated(val)) then
            call get_error(json, path, 'not an integer array')
        end if
        if (present(n)) call check_shape(path, shape(val), [n])
    end subroutine get_integer_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_real_1d
    !> @brief Get a real vector. With n, its size must be n.
    !----------------------------------------------------------------------------------------------
    subroutine get_real_1d(json, obj, path, val, n)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: val(:)
        integer, intent(in), optional :: n !< Expected size.
        type(json_value), pointer :: p

        p => json_get_child(json, obj, path)
        call vector_from(json, p, path, val)
        if (present(n)) call check_shape(path, shape(val), [n])
    end subroutine get_real_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_real_2d
    !> @brief Get a rank 2 real array written by add_real_2d.
    !> @details
    !! A nested array with no elements does not give the size of the inner dimensions. With
    !! expected, the shape must be expected, and an empty array gets this shape.
    !----------------------------------------------------------------------------------------------
    subroutine get_real_2d(json, obj, path, val, expected)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: val(:, :)
        integer, intent(in), optional :: expected(2) !< Expected shape.
        type(json_value), pointer :: p

        p => json_get_child(json, obj, path)
        call matrix_from(json, p, path, val)
        if (present(expected)) then
            if (size(val) == 0 .and. product(expected) == 0) then
                deallocate(val)
                allocate(val(expected(1), expected(2)))
            end if
            call check_shape(path, shape(val), expected)
        end if
    end subroutine get_real_2d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_real_3d
    !> @brief Get a rank 3 real array written by add_real_3d.
    !> @details
    !! With expected, the shape must be expected, and an empty array gets this shape.
    !----------------------------------------------------------------------------------------------
    subroutine get_real_3d(json, obj, path, val, expected)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: val(:, :, :)
        integer, intent(in), optional :: expected(3) !< Expected shape.
        type(json_value), pointer :: p, child
        real(dp), allocatable :: mat(:, :)
        integer :: i, n1

        p => json_get_child(json, obj, path)
        call check_array(json, p, path)
        n1 = json%count(p)
        do i = 1, n1
            call json%get_child(p, i, child)
            call matrix_from(json, child, path, mat)
            if (i == 1) allocate(val(n1, size(mat, 1), size(mat, 2)))
            call check_shape(path, shape(mat), shape(val(i, :, :)))
            val(i, :, :) = mat
        end do
        if (n1 == 0) allocate(val(0, 0, 0))
        if (present(expected)) then
            if (size(val) == 0 .and. product(expected) == 0) then
                deallocate(val)
                allocate(val(expected(1), expected(2), expected(3)))
            end if
            call check_shape(path, shape(val), expected)
        end if
    end subroutine get_real_3d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_complex_1d
    !> @brief Get a complex vector written by add_complex_1d. With n, its size must be n.
    !----------------------------------------------------------------------------------------------
    subroutine get_complex_1d(json, obj, path, val, n)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        complex(dp), allocatable, intent(out) :: val(:)
        integer, intent(in), optional :: n !< Expected size.
        type(json_value), pointer :: p
        real(dp), allocatable :: re(:), im(:)

        p => json_get_child(json, obj, path)
        call get_real_1d(json, p, 're', re)
        call get_real_1d(json, p, 'im', im, size(re))
        val = cmplx(re, im, kind=dp)
        if (present(n)) call check_shape(path, shape(val), [n])
    end subroutine get_complex_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_integer
    !> @brief Stop unless the integer at path equals expected, the value of the current run.
    !----------------------------------------------------------------------------------------------
    subroutine check_integer(json, obj, path, expected)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        integer, intent(in) :: expected
        integer :: val

        call get_integer(json, obj, path, val)
        if (val == expected) return
        write(stderr, '(1x,a)') 'Error in json_util_mod.'
        write(stderr, '(3x,3a,i0,a,i0,a)') 'Value "', path, '" is ', val, ', expected ', expected, &
        &                                  ' from the current setup.'
        stop 1
    end subroutine check_integer


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_string
    !> @brief Stop unless the string at path equals expected, the value of the current run.
    !----------------------------------------------------------------------------------------------
    subroutine check_string(json, obj, path, expected)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: expected
        character(len=:), allocatable :: val

        call get_string(json, obj, path, val)
        if (val == expected) return
        write(stderr, '(1x,a)') 'Error in json_util_mod.'
        write(stderr, '(3x,7a)') 'Value "', path, '" is "', val, '", expected "', expected,       &
        &                        '" from the current setup.'
        stop 1
    end subroutine check_string


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_logical
    !> @brief Stop unless the logical at path equals expected, the value of the current run.
    !----------------------------------------------------------------------------------------------
    subroutine check_logical(json, obj, path, expected)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj
        character(len=*), intent(in) :: path
        logical, intent(in) :: expected
        logical :: val

        call get_logical(json, obj, path, val)
        if (val .eqv. expected) return
        write(stderr, '(1x,a)') 'Error in json_util_mod.'
        write(stderr, '(3x,3a,l1,a,l1,a)') 'Value "', path, '" is ', val, ', expected ', expected, &
        &                                  ' from the current setup.'
        stop 1
    end subroutine check_logical


    !> Real vector from the array p.
    subroutine vector_from(json, p, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: p
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: val(:)

        call check_array(json, p, path)
        call json%get(p, val)
        if (json%failed() .or. .not. allocated(val)) call get_error(json, path, 'not a real array')
    end subroutine vector_from


    !> Rank 2 real array from the nested array p.
    subroutine matrix_from(json, p, path, val)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: p
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: val(:, :)
        type(json_value), pointer :: child
        real(dp), allocatable :: vec(:)
        integer :: i, n1

        call check_array(json, p, path)
        n1 = json%count(p)
        do i = 1, n1
            call json%get_child(p, i, child)
            call vector_from(json, child, path, vec)
            if (i == 1) allocate(val(n1, size(vec)))
            call check_shape(path, shape(vec), [size(val, 2)])
            val(i, :) = vec
        end do
        if (n1 == 0) allocate(val(0, 0))
    end subroutine matrix_from


    !> Stop unless p is an array.
    subroutine check_array(json, p, path)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: p
        character(len=*), intent(in) :: path
        integer :: var_type

        call json%info(p, var_type=var_type)
        if (var_type /= json_array) call get_error(json, path, 'not an array')
    end subroutine check_array


    !> Stop if the shape of the value at path is not the expected one.
    subroutine check_shape(path, found, expected)
        character(len=*), intent(in) :: path
        integer, intent(in) :: found(:)
        integer, intent(in) :: expected(:)

        if (size(found) == size(expected)) then
            if (all(found == expected)) return
        end if
        write(stderr, '(1x,a)') 'Error in json_util_mod.'
        write(stderr, '(3x,3a)') 'Value "', path, '" has the wrong shape.'
        write(stderr, '(3x,a,100(1x,i0))') 'Expected:', expected
        write(stderr, '(3x,a,100(1x,i0))') 'Found:   ', found
        stop 1
    end subroutine check_shape


    !> Stop with a message about the value at path.
    subroutine get_error(json, path, problem)
        type(json_core), intent(inout) :: json
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: problem
        character(len=:), allocatable :: msg
        logical :: ok

        write(stderr, '(1x,a)') 'Error in json_util_mod.'
        write(stderr, '(3x,5a)') 'Value "', path, '" ', problem, '.'
        call json%check_for_errors(ok, msg)
        if (.not. ok) write(stderr, '(3x,2a)') 'json-fortran: ', msg
        stop 1
    end subroutine get_error


end module json_util_mod
