!--------------------------------------------------------------------------------------------------
! MODULE: assignment_problem_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2017
!
!> @brief Subroutine for solving the assignment problem.
!> @details
!! The (linear) assignment problem:
!!
!! The problem consists of n agents and n tasks. Any agent can be assigned to perform any task,
!! incurring a cost depending on the agent-task assignment. It is required to perform all tasks
!! by assigning exactly one agent to each task and exactly one task to each agent in such a way
!! that the total cost (sum of agent-task costs) of the assignment is minimized.
!!
!! A modified form of the Hungarian method (Kuhn–Munkres algorithm) is used.
!! Adapted from: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
!!
!! References:
!! J. Munkres, "Algorithms for the Assignment and Transportation Problems", J. SIAM, 5(1),
!!              32–38, 1957 March.
!--------------------------------------------------------------------------------------------------
module assignment_problem_mod
    use global_defs
    implicit none

    private
    public :: assignment_problem
    public :: getloops
    public :: getloop


    real(dp), parameter :: tinydp = 1.e-15_dp

    integer :: m
    integer :: n
    integer :: k
    real(dp), allocatable :: cm(:, :)
    logical, allocatable :: zmat(:, :)
    logical, allocatable :: rowc(:)
    logical, allocatable :: colc(:)
    logical, allocatable :: star(:, :)
    logical, allocatable :: prim(:, :)
    integer :: cpos(2)
    integer :: step


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: assignment_problem
    !
    !> @brief The Kuhn-Munkres algorithm for square matrices.
    !> @details
    !! Input is a real matrix with dimensions n x n. The output is a vector containing the optimal
    !! permutation of columns associated to each row of the matrix (e. g. element 1 of the vector
    !! is the column associated with row 1 of the matrix).
    !!
    !! The algorithm is split into steps, each contained in a separate helper subroutine.
    !----------------------------------------------------------------------------------------------
    subroutine assignment_problem(matrix, assigned_rows, assigned_cols)
        real(dp), intent(in) :: matrix(:, :)
        integer, allocatable, intent(out) :: assigned_rows(:)
        integer, allocatable, intent(out) :: assigned_cols(:)

        if (size(matrix, 1) > size(matrix, 2)) then
            allocate(cm, source=transpose(matrix))
        else
            allocate(cm, source=matrix)
        end if
        m = size(cm, 1)
        n = size(cm, 2)
        k = min(m, n)
        allocate(zmat(m, n))
        allocate(star(m, n))
        allocate(prim(m, n))
        allocate(rowc(m))
        allocate(colc(n))
        if (.not. allocated(assigned_rows)) allocate(assigned_rows(k))
        if (.not. allocated(assigned_cols)) allocate(assigned_cols(k))

        step = 1
        colc = .false.
        rowc = .false.
        prim = .false.
        star = .false.
        do
            select case(step)
            case(1)
                call step1()
            case(2)
                call step2()
            case(3)
                call step3()
            case(4)
                call step4()
            case(5)
                call step5()
            case(6)
                call step6()
            case(7)
                if (size(matrix, 1) /= m) then
                    call step7(transpose(star), assigned_rows, assigned_cols)
                else
                    call step7(star, assigned_rows, assigned_cols)
                end if
                exit
            end select
        end do

        deallocate(cm)
        deallocate(zmat)
        deallocate(star)
        deallocate(prim)
        deallocate(rowc)
        deallocate(colc)
    end subroutine assignment_problem


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step1
    !
    !> @details
    !! For each row of the matrix, find the smallest element and subtract it from every element in
    !! its row. Repeat for columns. Create a logical matrix which is .true. where the value of the
    !! input matrix is zero. Go to Step 2.
    !----------------------------------------------------------------------------------------------
    subroutine step1()
        integer :: i
        do i = 1, m
            cm(i, :) = cm(i, :) - minval(cm(i, :))
        end do
        zmat = (cm < tinydp)
        step = 2
    end subroutine step1


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step2
    !
    !> @details
    !! Cycle through the zeroes. If there is no starred zero in the row or column of a zero, star
    !! it. Go to Step 3.
    !----------------------------------------------------------------------------------------------
    subroutine step2()
        integer :: i, j
        star = .false.
        do i = 1, m
            if (any(star(i, :))) cycle
            do j = 1, n
                if (.not. zmat(i, j)) cycle
                if (any(star(:, j))) cycle
                star(i, j) = .true.
                exit
            end do
        end do
        step = 3
    end subroutine step2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step3
    !
    !> @details
    !! Cover each column containing a starred zero. If k columns are covered, the starred zeros
    !! describe a complete set of unique assignments (go to step 7). If not, go to step 4.
    !----------------------------------------------------------------------------------------------
    subroutine step3()
        integer :: j, c

        c = 0
        do j = 1, n
            if (.not. any(star(:, j))) cycle
            colc(j) = .true.
            c = c + 1
        end do

        if (c == k) then
            step = 7
        else
            step = 4
        end if
    end subroutine step3


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step4
    !
    !> @details
    !! Find a noncovered zero and prime it.  If there is no starred zero in the row containing the
    !! primed zero, mark the position of the primed zero and Go to Step 5.  Otherwise, cover this
    !! row and uncover the column containing the starred zero. Repeat. If there are no uncovered
    !! zeros left, go to Step 6.
    !----------------------------------------------------------------------------------------------
    subroutine step4()
        integer :: i, j
        logical :: starincol

        do
            cpos = [0, 0]
            starincol = .false.
outer:      do i = 1, m
                if (rowc(i)) cycle
inner:          do j = 1, n
                    if (colc(j)) cycle
                    if (.not. zmat(i, j)) cycle
                    prim(i, j) = .true.
                    cpos = [i, j]
                    exit outer
                end do inner
            end do outer
            if (cpos(1) == 0) then
                step = 6
                return
            end if
            do j = 1, n
                if (.not. star(cpos(1), j)) cycle
                colc(j) = .false.
                rowc(cpos(1)) = .true.
                starincol = .true.
            end do
            if (.not. starincol) then
                step = 5
                return
            end if
        end do

    end subroutine step4


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step5
    !
    !> @details
    !! Starting at the final primed zero of Step 4, move through the matrix in the following way:
    !! Go to a starred zero in the same column as the primed zero (if possible). Go to a primed
    !! zero in the same row as the starred zero (always possible). Along the way, unstar all
    !! starred zeros and star all primed zeros. When no starred zeros are present in the column of
    !! a primed zero, uncover all lines of the matrix, unprime all zeros and go to Step 3.
    !----------------------------------------------------------------------------------------------
    subroutine step5()
        integer :: i, j
        logical :: starincol

        do
            starincol = .false.
            do i = 1, m
                if (.not. star(i, cpos(2))) cycle
                starincol = .true.
                star(cpos(1), cpos(2)) = .true.
                cpos(1) = i
                exit
            end do
            if (.not. starincol) then
                star(cpos(1), cpos(2)) = .true.
                prim = .false.
                rowc = .false.
                colc = .false.
                step = 3
                return
            end if
            do j = 1, n
                if (.not. prim(cpos(1), j)) cycle
                star(cpos(1), cpos(2)) = .false.
                cpos(2) = j
                exit
            end do
        end do
    end subroutine step5


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step6
    !
    !> @details
    !! Find the smallest uncovered value. Add this value to every element of each covered row, and
    !! subtract it from every element of each uncovered column. Return to Step 4.
    !----------------------------------------------------------------------------------------------
    subroutine step6()
        real(dp) :: mval
        logical :: mask(m, n)
        integer :: i, j

        mask = .true.
        do i = 1, m
            if (rowc(i)) mask(i, :) = .false.
        end do
        do j = 1, n
            if (colc(j)) mask(:, j) = .false.
        end do
        mval = minval(cm, mask)
        do i = 1, m
            if (rowc(i)) cm(i, :) = cm(i, :) + mval
        end do
        do j = 1, n
            if (.not. colc(j)) cm(:, j) = cm(:, j) - mval
        end do
        zmat = (cm < tinydp)
        step = 4
    end subroutine step6


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step7
    !
    !> @details
    !! Construct the output vector from the matrix of starred zeros.
    !----------------------------------------------------------------------------------------------
    subroutine step7(solution, assigned_rows, assigned_cols)
        logical, intent(in) :: solution(:, :)
        integer, intent(out) :: assigned_rows(:)
        integer, intent(out) :: assigned_cols(:)
        integer :: i, j, c

        c = 0
        do i = 1, size(solution, 1)
            do j = 1, size(solution, 2)
                if (solution(i, j)) then
                    c = c + 1
                    assigned_rows(c) = i
                    assigned_cols(c) = j
                end if
            end do
        end do
    end subroutine step7


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: GetLoops
    !> @brief Find all loops in given list of indexes.
    !----------------------------------------------------------------------------------------------
    subroutine getloops(list, loops)
        use ivec_mod
        integer, intent(in) :: list(:)
        type(ivec), allocatable, intent(out) :: loops(:)
        type(ivec) :: tmp(size(list))
        logical :: mask(size(list))
        integer :: i, n

        mask = .true.
        n = 0
        do i = 1, size(list)
            if (.not. mask(i)) cycle
            n = n + 1
            call getloop(list, i, tmp(n)%c, mask)
        end do

        if (allocated(loops)) deallocate(loops)
        allocate(loops(n))
        do i = 1, n
            loops(i) = tmp(i)
        end do
    end subroutine getloops


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: GetLoop
    !> @brief Follow list of indexes until it loops to initial index.
    !----------------------------------------------------------------------------------------------
    subroutine getloop(list, init, loop, mask)
        integer, intent(in) :: list(:) !< List of indexes to follow.
        integer, intent(in) :: init !< Initial position in list.
        integer, allocatable, intent(out) :: loop(:) !< List of indexes until init is found.
        logical :: mask(size(list))
        integer :: i, n

        n = 1
        i = init
        getn: do
            mask(i) = .false.
            if (list(i) == init) exit getn
            i = list(i)
            n = n + 1
            if (i > size(list)) then
                write(stderr, *) 'Error in getloop subroutine. Index out of bounds.'
                call abort()
            end if
            if (n > size(list)) then
                write(stderr, *) 'Error in getloop subroutine. Bad input list.'
                call abort()
            end if
        end do getn

        if (allocated(loop)) deallocate(loop)
        allocate(loop(n))
        i = init
        fill: do n = 1, size(loop)
            loop(n) = i
            if (list(i) == init) exit fill
            i = list(i)
        end do fill
    end subroutine getloop


end module assignment_problem_mod
