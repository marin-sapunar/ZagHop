program util_test
    implicit none
    integer :: ret


    ret = assignment_test()
    if (ret /= 0) then
        write(*,*) "Error in assignment_test."
        stop ret
    end if

    ret = random_test()
    if (ret /= 0) then
        write(*,*) "Error in random_test."
        stop ret
    end if

    ret = string_test()
    if (ret /= 0) then
        write(*,*) "Error in string_test."
        stop ret
    end if

    contains
    
    
    integer function random_test() result(ret)
        use global_defs, only : dp
        use random_mod, only : rng_type
        use random_pcg_mod, only : rng_pcg_xsh_rr
        class(rng_type), allocatable :: rng
        real(dp) :: rnum(1000)
        integer :: i
        real(dp) :: tiny = 1.0e-15

        ret = 1
        allocate(rng_pcg_xsh_rr::rng)
        call rng%init(100)
        do i = 1, 100
            call rng%uniform(rnum)
        end do
        if (abs(rnum(1)-0.84574914766609421_dp) > tiny) return
        if (abs(rnum(999)-0.32289557934806112_dp) > tiny) return
        ret = 0
    end function random_test
 

    integer function assignment_test() result(ret)
        use global_defs, only : dp
        use ivec_mod, only : ivec
        use assignment_problem_mod, only : assignment_problem, getloops
        real(dp) :: matrix(4, 3)
        integer, allocatable :: rows(:)
        integer, allocatable :: cols(:)
        type(ivec), allocatable :: loops(:)

        matrix(1, :) = [8.0_dp, 4.0_dp, 7.0_dp]
        matrix(2, :) = [5.0_dp, 2.0_dp, 3.0_dp]
        matrix(3, :) = [9.0_dp, 6.0_dp, 7.0_dp]
        matrix(4, :) = [9.0_dp, 4.0_dp, 8.0_dp]

        ret = 1
        call assignment_problem(matrix, rows, cols)
        if (any(rows /= [1, 2, 4])) return
        if (any(cols /= [1, 3, 2])) return
        call assignment_problem(transpose(matrix), rows, cols)
        if (any(rows /= [1, 2, 3])) return
        if (any(cols /= [1, 4, 2])) return

        matrix(1, :) = [0.0_dp, 0.0_dp, 1.0_dp]
        matrix(2, :) = [0.0_dp, 1.0_dp, 0.0_dp]
        matrix(3, :) = [1.0_dp, 0.0_dp, 0.0_dp]
        matrix(4, :) = [0.0_dp, 0.0_dp, 0.0_dp]
        call assignment_problem(-matrix(1:3, :), rows, cols)
        if (any(rows /= [1, 2, 3])) return
        if (any(cols /= [3, 2, 1])) return
        call getloops(cols, loops)
        if (size(loops) /= 2) return
        if (any(loops(1)%c /= [1, 3])) return
        if (any(loops(2)%c /= [2])) return
        ret = 0
    end function assignment_test


    integer function string_test() result(ret)
        use string_mod
        integer :: narg
        type(string), allocatable :: args(:)
        integer, allocatable :: indexlist(:)

        ret = 1
        call string_parse("a=5", "= ", narg, args)
        if (narg /= 2) return
        if (args(1)%s /= "a") return
        call string_parse("a = 5, b 7, d+= 500", "=, ", narg, args)
        if (narg /= 6) return
        if (args(1)%s /= "a") return
        if (args(6)%s /= "500") return
        call read_index_list("1,2-5, 300,1000-5000", indexlist)
        if (size(indexlist) /= 4007) return
        if (indexlist(6) /= 300) return
        if (indexlist(4006) /= 4999) return
        if (tolower("c") /= "c") return
        if (tolower("C") /= "c") return
        if (tolower("Ca") /= "ca") return
        if (toupper("Ca") /= "CA") return
        ret = 0
    end function string_test

    
end program util_test
