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

    ret = rigid_body_test()
    if (ret /= 0) then
        write(*,*) "Error in rigid_body_test."
        stop ret
    end if

    ret = kabsch_derivative_test()
    if (ret /= 0) then
        write(*,*) "Error in kabsch_derivative_test."
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


    integer function rigid_body_test() result(ret)
        use global_defs, only : dp
        use matrix_mod, only : unit_mat
        use rigid_body_mod, only : external_modes, kabsch_rotation
        real(dp) :: mass(4)
        real(dp) :: geom(3, 4)
        real(dp) :: moved(3, 4)
        real(dp) :: lin(3, 3)
        real(dp) :: axis(3), cross(3, 3), rot(3, 3), u(3, 3)
        real(dp), allocatable :: ext(:, :)
        real(dp) :: tiny = 1.0e-12_dp

        mass = [12.0_dp, 1.0_dp, 16.0_dp, 14.0_dp]
        geom = reshape([0.1_dp, 0.2_dp, -0.3_dp, 2.1_dp, 0.1_dp, 0.4_dp, -0.4_dp, 2.2_dp,        &
        &               0.3_dp, 0.5_dp, -1.1_dp, 1.9_dp], [3, 4])
        ! Rotation by 2 rad about the (1, 2, 2) axis.
        axis = [1.0_dp, 2.0_dp, 2.0_dp] / 3.0_dp
        cross = reshape([0.0_dp, axis(3), -axis(2), -axis(3), 0.0_dp, axis(1), axis(2),          &
        &                -axis(1), 0.0_dp], [3, 3])
        rot = cos(2.0_dp) * unit_mat(3) + sin(2.0_dp) * cross                                    &
        &   + (1.0_dp - cos(2.0_dp)) * spread(axis, 2, 3) * spread(axis, 1, 3)

        ret = 1
        ! Orthonormal translation/rotation vectors: 6 for a nonlinear molecule, 3 without
        ! rotation and 5 for a linear molecule.
        ext = external_modes(mass, geom, .true.)
        if (size(ext, 2) /= 6) return
        if (maxval(abs(matmul(transpose(ext), ext) - unit_mat(6))) > tiny) return
        ext = external_modes(mass, geom, .false.)
        if (size(ext, 2) /= 3) return
        lin = 0.0_dp
        lin(3, :) = [-2.2_dp, 0.0_dp, 2.2_dp]
        ext = external_modes(mass(1:3), lin, .true.)
        if (size(ext, 2) /= 5) return
        if (maxval(abs(matmul(transpose(ext), ext) - unit_mat(5))) > tiny) return

        ! Kabsch rotation undoes the rotation of a translated geometry.
        moved = matmul(rot, geom) + spread([1.0_dp, -2.0_dp, 0.5_dp], 2, 4)
        u = kabsch_rotation(mass, moved, geom)
        if (maxval(abs(matmul(u, rot) - unit_mat(3))) > tiny) return
        ! Planar molecule (singular covariance matrix), the result should still be a rotation.
        geom(1, :) = 0.0_dp
        moved = matmul(rot, geom)
        u = kabsch_rotation(mass, moved, geom)
        if (maxval(abs(matmul(u, rot) - unit_mat(3))) > tiny) return
        ret = 0
    end function rigid_body_test


    !> Compare kabsch_derivative with finite differences of f(rho(R)), where rho are the aligned
    !! positions and f(rho) = sum_i c_i . (rho_i - r0_i) + |rho - r0|^2 / 2.
    integer function kabsch_derivative_test() result(ret)
        use global_defs, only : dp
        use rigid_body_mod, only : center_of_mass, kabsch_rotation, kabsch_derivative
        real(dp) :: mass(4)
        real(dp) :: ref(3, 4)
        real(dp) :: geom(3, 4)
        real(dp) :: coef(3, 4)
        real(dp) :: deriv(3, 4), num(3, 4), r(3, 4)
        real(dp) :: u(3, 3)
        real(dp) :: torque(3)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp), parameter :: tol = 1.0e-8_dp
        real(dp), parameter :: tiny = 1.0e-12_dp
        integer :: i, k

        mass = [12.0_dp, 1.0_dp, 16.0_dp, 14.0_dp]
        ref = reshape([0.1_dp, 0.2_dp, -0.3_dp, 2.1_dp, 0.1_dp, 0.4_dp, -0.4_dp, 2.2_dp,         &
        &              0.3_dp, 0.5_dp, -1.1_dp, 1.9_dp], [3, 4])
        coef = reshape([0.3_dp, -0.2_dp, 0.5_dp, 0.1_dp, 0.4_dp, -0.6_dp, -0.3_dp, 0.2_dp,       &
        &               0.1_dp, 0.2_dp, -0.1_dp, 0.3_dp], [3, 4])
        ! Displaced geometry, rotated by 90 degrees about z and translated.
        geom(1, :) = ref(2, :) + 0.2_dp * coef(1, :) + 1.0_dp
        geom(2, :) = -ref(1, :) + 0.2_dp * coef(2, :) - 0.5_dp
        geom(3, :) = ref(3, :) + 0.2_dp * coef(3, :)

        ret = 1
        u = kabsch_rotation(mass, geom, ref)
        r = matmul(u, geom - spread(center_of_mass(mass, geom), 2, 4))                           &
        &   - (ref - spread(center_of_mass(mass, ref), 2, 4))
        deriv = kabsch_derivative(mass, geom, ref, u, coef + r)
        do i = 1, 4
            do k = 1, 3
                r = geom
                r(k, i) = r(k, i) + step
                num(k, i) = aligned_test_function(mass, ref, coef, r)
                r(k, i) = r(k, i) - 2 * step
                num(k, i) = num(k, i) - aligned_test_function(mass, ref, coef, r)
            end do
        end do
        num = num / (2 * step)
        if (maxval(abs(deriv - num)) > tol) return
        ! No net force or torque.
        if (maxval(abs(sum(deriv, dim=2))) > tiny) return
        r = geom - spread(center_of_mass(mass, geom), 2, 4)
        torque = 0.0_dp
        do i = 1, 4
            torque = torque + [r(2, i) * deriv(3, i) - r(3, i) * deriv(2, i),                    &
            &                  r(3, i) * deriv(1, i) - r(1, i) * deriv(3, i),                    &
            &                  r(1, i) * deriv(2, i) - r(2, i) * deriv(1, i)]
        end do
        if (maxval(abs(torque)) > tiny) return
        ret = 0
    end function kabsch_derivative_test


    function aligned_test_function(mass, ref, coef, geom) result(val)
        use global_defs, only : dp
        use rigid_body_mod, only : center_of_mass, kabsch_rotation
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: ref(:, :)
        real(dp), intent(in) :: coef(:, :)
        real(dp), intent(in) :: geom(:, :)
        real(dp) :: val
        real(dp) :: disp(3, size(mass))

        disp = matmul(kabsch_rotation(mass, geom, ref),                                           &
        &             geom - spread(center_of_mass(mass, geom), 2, size(mass)))                  &
        &    - (ref - spread(center_of_mass(mass, ref), 2, size(mass)))
        val = sum(coef * disp) + 0.5_dp * sum(disp**2)
    end function aligned_test_function

    
end program util_test
