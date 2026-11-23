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

    ret = rng_state_test()
    if (ret /= 0) then
        write(*,*) "Error in rng_state_test."
        stop ret
    end if

    ret = json_util_test()
    if (ret /= 0) then
        write(*,*) "Error in json_util_test."
        stop ret
    end if

    ret = truncate_file_test()
    if (ret /= 0) then
        write(*,*) "Error in truncate_file_test."
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


    !> After loading a saved state, a fresh generator continues the sequence of the original one,
    !! including the cached second number of a Gaussian pair.
    integer function rng_state_test() result(ret)
        use global_defs, only : dp
        use json_module, only : json_core, json_value
        use random_mod, only : rng_type
        use random_pcg_mod, only : rng_pcg_xsh_rr
        class(rng_type), allocatable :: rng, copy
        type(json_core) :: json
        type(json_value), pointer :: root, back
        character(len=:), allocatable :: text
        real(dp) :: r1(50), r2(50), g1(3), g2(3)

        ret = 1
        allocate(rng_pcg_xsh_rr::rng)
        allocate(rng_pcg_xsh_rr::copy)
        call rng%init(4242)
        call rng%uniform(r1)
        call rng%gaussian(0.0_dp, 1.0_dp, g1(1))
        call json%initialize()
        call json%create_object(root, '')
        call rng%save_state(json, root)
        call json%serialize(root, text)
        call json%destroy(root)
        call rng%gaussian(0.0_dp, 1.0_dp, g1(2:3))
        call rng%uniform(r1)

        call copy%init(1)
        call json%deserialize(back, text)
        call copy%load_state(json, back)
        call json%destroy(back)
        call copy%gaussian(0.0_dp, 1.0_dp, g2(2:3))
        call copy%uniform(r2)
        if (copy%seed /= 4242) return
        if (any(g1(2:3) /= g2(2:3))) return
        if (any(r1 /= r2)) return
        ret = 0
    end function rng_state_test


    !> Arrays written with json_util_mod and read back from the serialized text are unchanged.
    integer function json_util_test() result(ret)
        use global_defs, only : dp
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array, json_get_required, json_get_child
        type(json_core) :: json
        type(json_value), pointer :: root, back, obj
        character(len=:), allocatable :: text, str
        real(dp) :: a2(3, 2), a3(2, 3, 4), x
        complex(dp) :: c1(3)
        real(dp), allocatable :: b2(:, :), b3(:, :, :), e2(:, :), v(:)
        complex(dp), allocatable :: d1(:)
        integer, allocatable :: iv(:)
        integer :: i, j, k, n

        ret = 1
        do i = 1, 3
            do j = 1, 2
                a2(i, j) = 1.0_dp / (i + 7 * j) - 0.1_dp * i
            end do
        end do
        do i = 1, 2
            do j = 1, 3
                do k = 1, 4
                    a3(i, j, k) = sqrt(real(i + 10 * j + 100 * k, dp)) * 1.0e-7_dp
                end do
            end do
        end do
        c1 = cmplx([1.0_dp / 3.0_dp, -2.0_dp, 0.0_dp], [0.1_dp, 0.0_dp, -1.0e-300_dp], kind=dp)

        call json%initialize()
        call json%create_object(root, '')
        call json%create_object(obj, 'section')
        call json%add(root, obj)
        call json_add_array(json, obj, 'a2', a2)
        call json_add_array(json, obj, 'a3', a3)
        call json_add_array(json, obj, 'c1', c1)
        call json_add_array(json, obj, 'empty', reshape([real(dp) ::], [0, 4]))
        call json%add(obj, 'x', 4.0_dp * atan(1.0_dp))
        call json%add(obj, 'n', -7)
        call json%add(obj, 'iv', [3, 1, 2])
        call json%add(obj, 'str', 'text')
        call json%serialize(root, text)
        call json%destroy(root)

        call json%deserialize(back, text)
        if (json%failed()) return
        obj => json_get_child(json, back, 'section')
        call json_get_required(json, obj, 'a2', b2)
        if (any(shape(b2) /= shape(a2))) return
        if (any(b2 /= a2)) return
        call json_get_required(json, obj, 'a3', b3, shape(a3))
        if (any(b3 /= a3)) return
        call json_get_required(json, obj, 'c1', d1, 3)
        if (any(d1 /= c1)) return
        call json_get_required(json, obj, 'empty', e2, [0, 4])
        if (any(shape(e2) /= [0, 4])) return
        call json_get_required(json, obj, 'x', x)
        if (x /= 4.0_dp * atan(1.0_dp)) return
        call json_get_required(json, obj, 'n', n)
        if (n /= -7) return
        call json_get_required(json, obj, 'iv', iv, 3)
        if (any(iv /= [3, 1, 2])) return
        call json_get_required(json, obj, 'str', str)
        if (str /= 'text') return
        call json_get_required(json, obj, 'a2(2)', v, 2)
        if (any(v /= a2(2, :))) return
        call json%destroy(back)
        ret = 0
    end function json_util_test


    !> A file is cut to the given number of bytes and can be appended to afterwards.
    integer function truncate_file_test() result(ret)
        use global_defs, only : int64
        use file_mod, only : truncate_file
        character(len=*), parameter :: fname = 'util_test_truncate.txt'
        character(len=20) :: line
        integer(int64) :: n_byte
        integer :: u, i, ios

        ret = 1
        open(newunit=u, file=fname, action='write', status='replace')
        write(u, '(a)') 'line 1'
        flush(u)
        inquire(unit=u, size=n_byte)
        do i = 2, 5
            write(u, '(a,i0)') 'line ', i
        end do
        close(u)
        call truncate_file(fname, n_byte)
        open(newunit=u, file=fname, action='write', position='append')
        write(u, '(a)') 'line 6'
        close(u)
        open(newunit=u, file=fname, action='read')
        read(u, '(a)') line
        if (line /= 'line 1') return
        read(u, '(a)') line
        if (line /= 'line 6') return
        read(u, '(a)', iostat=ios) line
        if (ios == 0) return
        close(u, status='delete')
        ret = 0
    end function truncate_file_test

    
end program util_test
