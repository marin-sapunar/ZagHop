!--------------------------------------------------------------------------------------------------
! MODULE: interface_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date January 2020
!--------------------------------------------------------------------------------------------------
module interface_mod
    use global_defs
    use control_var
    implicit none

    private
    public :: run_qm
    public :: run_mm


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Interface_QMRun
    !
    ! DESCRIPTION:
    !> @brief Run a single point QM calculation.
    !> @details
    !! @todo UPDATE DOCS
    !----------------------------------------------------------------------------------------------
    subroutine run_qm(t, hop)
#ifdef QUANTICS
        use shzagreb_inter, only : shzagreb_run
#endif
        use system_type_mod, only : system_type
        use matrix_mod, only : unit_mat
        use model_mod, only : qmodel
        use json_module, only : json_core, json_file, json_value
        type(system_type), intent(inout) :: t
        logical, intent(in) :: hop
        integer :: cunit, i, j, d1, d2, d3, d4
        logical :: check(5)
        character(len=200) :: json_str
        real(dp), allocatable :: rn1(:)
        real(dp), allocatable :: rn2(:, :)
        real(dp), allocatable :: rn3(:, :, :)
        type(json_core) :: json
        type(json_file) :: read_json
        type(json_value), pointer :: base, p_top, p_array_i, p_array_j
        type(json_value), pointer :: p_array_0, p_geom_entry
        real(dp), allocatable :: wrk_en(:)
        real(dp), allocatable :: wrk_nadv(:, :, :)
        real(dp), allocatable :: wrk_soc(:, :)
        real(dp), allocatable :: wrk_adt(:, :)
        integer, allocatable :: spinv(:)
        logical :: found

        !> @todo Move this allocation check to wf object?
        do i = 1, t%wf%n_state
            if (t%wf%need_gradient(i)) then
                if (.not. allocated(t%wf%qm_state(i)%gradient)) then
                    allocate(t%wf%qm_state(i)%gradient(t%ndim, t%qnatom))
                end if
            end if
            do j = 1, t%wf%n_state
                if (t%wf%need_nadv(i, j)) then
                    if (.not. allocated(t%wf%qm_state(i)%nadv(j)%c)) then
                        allocate(t%wf%qm_state(i)%nadv(j)%c(t%ndim * t%qnatom))
                    end if
                end if
            end do
            if (t%wf%need_soc(i)) then
                if (.not. allocated(t%wf%qm_state(i)%soc)) then
                    allocate(t%wf%qm_state(i)%soc(t%wf%n_state))
                end if
            end if
        end do

        select case (ctrl%qlib)
        case(0)
            if ((t%step > 0) .and. (.not. hop)) then
                call system('rm -rf prevstep')
                call system('cp -r '//ctrl%qmdir//' prevstep')
            end if

            call json%initialize()
            call json%create_object(base, "")

            call json%create_object(p_top, "step")
            call json%add(base, p_top)
            call json%add(p_top, "step", t%step)

            call json%create_array(p_top, "system")
            call json%add(base, p_top)
            !@todo Loop over subsystems for QM/MM.
            call json%create_object(p_array_i, "")
            call json%add(p_top, p_array_i)
            call json%add(p_array_i, "natom", t%qnatom)
            call json%create_array(p_array_j, "geom")
            call json%add(p_array_i, p_array_j)
            do j = 1, t%qnatom
                call json%add(p_array_j, "", t%geom(1, t%qind(j)))
                call json%add(p_array_j, "", t%geom(2, t%qind(j)))
                call json%add(p_array_j, "", t%geom(3, t%qind(j)))
            end do
            ! done

            call json%create_array(p_top, "states")
            call json%add(base, p_top)
            do i = 1, t%wf%n_state_group
                ! Create entry for each state group.
                call json%create_object(p_array_i, "")
                call json%add(p_top, p_array_i)
                call json%add(p_array_i, "system", 1)
                call json%add(p_array_i, "nstate", t%wf%n_state_per_group(i))
                call json%add(p_array_i, "multiplicity", t%wf%qm_state(t%wf%index(i, 1))%multiplicity)
                call json%add(p_array_i, "energy", .true.)
                call json%add(p_array_i, "oscillator_strength", ctrl%oscill)
                call json%create_array(p_array_j, "gradient")
                call json%add(p_array_i, p_array_j)
                do j = 1, t%wf%n_state_per_group(i)
                    if (t%wf%index(i, j) == t%wf%active_state) call json%add(p_array_j, "", j)
                end do
                !> @todo Add requests for NACs and SOCs if needed.
            end do
            call json%print(base, "qm.json")
            call json%destroy(base)

            if (ctrl%mm) then
                open(newunit=cunit, file='mm_geom', action='write')
                do i = 1, t%mnatom
                    write(cunit, *) t%geom(:, t%mind(i))
                end do
                close(cunit)
            end if

            ! Call interface.
            call system('rm -f qm_out.json')
            call system(ctrl%qprog)

            ! Check if energy and gradient files were created.
            inquire(file='qm_out.json', exist=check(1))
            if (.not. check(1)) then
                call errstop("run_qm", "File qm_out.json not found after QM calculation.", 1)
            end if

            call read_json%load_file("qm_out.json")
            do i = 1, t%wf%n_state
                write(json_str, '(a,i0,a,i0,a)') "states(", t%wf%qm_state(i)%group, ").energy(", &
                &                                t%wf%qm_state(i)%group_state, ")"

                call read_json%get(json_str, t%wf%qm_state(i)%energy, found)
                if (.not. found) then
                    call errstop("run_qm", "Energy not found in QM output.", i)
                end if
            end do
            call json_get_2d(read_json, "states(1).gradient(1)", t%wf%qm_state(t%wf%active_state)%gradient)
            if (ctrl%oscill) then
            !    call read_json%get("states(1).oscillator_strength", t%qo)
                if (.not. found) then
                    call errstop("run_qm", "Oscillator strength not found in QM output.", 1)
                end if
            end if

            ! @todo these loops are not functional for now since the files can have only
            ! one set of couplings. Need to add this information to the json input and output.
            ! do i = 1, t%wf%n_state_group
            !     do j = 1, t%wf%n_state_group
            !         if ((t%couplings(i, j) == "overlap") .and. t%step /= 0) then
            !             call system(ctrl%oprog)
            !             open(newunit=cunit, file='qm_olap', action='read')
            !             read(cunit, *) d1, d2
            !             allocate(t%wf%overlap(1)%c(d1, d2))
            !             do i = 1, d1
            !                 read(cunit, *) t%wf%overlap(1)%c(i, :)
            !             end do
            !             close(cunit)
            !             ! @todo Check dimensions for variable_nstate runs, previous versions
            !             ! had min(d1, d2) for the loop and set the rest to a unit matrix.
            !         else if (t%couplings(i, j) == "nadvec") then
            !             open(newunit=cunit, file='qm_nadvec', action='read')
            !             do i = 1, t%max_nstate
            !                 do j = 1, t%max_nstate
            !                     read(cunit, *) t%couplings(i, j)%term(:, i, j)
            !                 end do
            !             end do
            !             close(cunit)
            !         end if
            !     end do
            ! end do
        case(1)
#ifdef QUANTICS
            !> @todo This is a workaround which avoids changing quantics_inter.f90 for now,
            ! but the interface should be modified so temporary arrays are not required.
            d1 = t%wf%n_state
            allocate(wrk_en(d1))
            allocate(wrk_nadv(t%natom * t%ndim, d1, d1))
            allocate(wrk_soc(d1, d1))
            spinv = t%wf%qm_state(:)%multiplicity
            call shzagreb_run(t%step, t%geom, t%wf%active_state, wrk_en, t%wf%qm_state(t%wf%active_state)%gradient, wrk_nadv, &
            &                 wrk_soc, spinv, ctrl%socbas, t%adt)
            t%wf%qm_state(:)%energy = wrk_en
            do i = 1, t%wf%n_state
                do j = 1, t%wf%n_state
                    if (i == j) cycle
                    if (t%wf%qm_state(i)%group == t%wf%qm_state(j)%group) then
                        d1 = d1 + 1
                        t%wf%qm_state(i)%nadv(j)%c = wrk_nadv(:, i, j)
                    end if
                    if (t%wf%qm_state(i)%multiplicity /= t%wf%qm_state(j)%multiplicity) then
                        t%wf%qm_state(i)%soc(j) = wrk_soc(i, j)
                    end if
                end do
            end do
#else
            call errstop("run_qm", "Code not compiled with quantics interface.", 1)
#endif
        case(2)
            !> @todo This is a workaround which avoids changing qmodel%eval for now,
            ! but the interface should be modified so temporary arrays are not required.
            wrk_en = t%wf%qm_state(:)%energy
            if (ctrl%tdc_type == 2 .or. ctrl%vrescale == 3) then
                allocate(wrk_nadv(t%natom * t%ndim, t%wf%n_state, t%wf%n_state))
            end if
            if (ctrl%tdc_type == 3) then
                allocate(wrk_adt(t%wf%n_state, t%wf%n_state))
            end if
            if (.not. allocated(t%wf%qm_state(t%wf%active_state)%gradient)) then
                allocate(t%wf%qm_state(t%wf%active_state)%gradient(t%ndim, t%qnatom))
            end if
            call qmodel%eval(t%geom, t%wf%active_state, wrk_en, t%wf%qm_state(t%wf%active_state)%gradient, &
            &                wrk_nadv, wrk_adt)
            t%wf%qm_state(:)%energy = wrk_en
            if ((ctrl%tdc_type == 2 ) .or. (ctrl%vrescale == 3)) then
                do i = 1, t%wf%n_state
                    do j = 1, t%wf%n_state
                        if (.not. t%wf%need_nadv(i, j)) cycle
                        t%wf%qm_state(i)%nadv(j)%c = wrk_nadv(:, i, j)
                    end do
                end do
            end if
            if (ctrl%tdc_type == 3) then
                t%wf%overlap(2)%c = wrk_adt
            end if
        case default
            call errstop("run_qm", "Unrecognized QM interface.", 1)
        end select

        do i = 1, t%qnatom
            t%grad(:, t%qind(i)) = t%wf%qm_state(t%wf%active_state)%gradient(:, i)
        end do

        ! Add random noise to evaluated values if requested.
        if (ctrl%noise > 0.0_dp) then
            allocate(wrk_en(t%wf%n_state))
            call ctrl%rng%uniform(wrk_en)
            wrk_en = (wrk_en - 0.5_dp) * ctrl%noise
            t%wf%qm_state(:)%energy = t%wf%qm_state(:)%energy + wrk_en
        end if
    end subroutine run_qm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Interface_MMRun
    !----------------------------------------------------------------------------------------------
    subroutine run_mm(t)
        use system_type_mod, only : system_type
        type(system_type), intent(inout) :: t
        integer :: i, cunit
        logical :: check1, check2, check3

        call nad_interface_write_atoms(t%natom, t%qnatom, t%mnatom, t%qind, t%mind, t%sym, t%geom,&
        &                              t%velo, ctrl%pbc, t%pbcbox)

        ! Call QM program.
        call system('rm -f mm_energy mm_grad mm_pbc')
        call system(ctrl%mprog)

        ! Check if energy and gradient files were created.
        inquire(file='mm_energy', exist=check1)
        if (.not. check1) write(stderr, *) 'Error, mm_energy file not found after MM calculation.'
        inquire(file='mm_grad', exist=check2)
        if (.not. check2) write(stderr, *) 'Error, mm_grad file not found after MM calculation.'
        check3 = .true.
        if (ctrl%pbc) then
            inquire(file='mm_pbc', exist=check3)
            if (.not. check3) write(stderr, *) 'Error, mm_pbc file not found after MM calculation.'
        end if
        if ((.not. check1) .or. (.not. check2) .or. (.not. check3)) then
            call system('rm -rf mmdir_error')
            call system('cp -r '//ctrl%qmdir//' mmdir_error')
            stop
        end if

        open(newunit=cunit, file='mm_energy', action='read')
        read(cunit, *) t%me
        close(cunit)

        open(newunit=cunit, file='mm_grad', action='read')
        do i = 1, t%natom
            read(cunit, *) t%grad(:, i)
        end do
        close(cunit)

        open(newunit=cunit, file='mm_pbc', action='read')
        read(cunit, *) t%pbcbox
        close(cunit)

        t%me(1) = t%me(1) - eelp(t%geom(:, t%qind), t%geom(:, t%mind), t%chrg(t%qind),             &
        &         t%chrg(t%mind), ctrl%mmcut)
    end subroutine run_mm


    subroutine nad_interface_write_atoms(nat, qnat, mnat, qind, mind, sym, geo, vel, pbc, box)
        integer, intent(in) :: nat, qnat, mnat
        integer, intent(in) :: qind(:), mind(:)
        character(len=2), intent(in) :: sym(:)
        real(dp), intent(in) :: geo(:, :), vel(:, :)
        logical, intent(in) :: pbc
        real(dp), intent(in) :: box(6)
        character(len=*), parameter :: ifmt = '(1000(i0,1x))'
        character(len=*), parameter :: atomfmt = '(1x,a2,2x,1000e24.16)'
        character(len=*), parameter :: relfmt = '(1000(e24.16,2x))'
        integer :: ounit, i

        open(newunit=ounit, file='mm_data', action='write')
        write(ounit, *) nat, qnat, mnat
        write(ounit, ifmt) qind
        write(ounit, ifmt) mind
        do i = 1, nat
            write(ounit, atomfmt) sym(i), geo(:, i), vel(:, i)
        end do
        if (pbc) write(ounit, relfmt) box
        close(ounit)
    end subroutine nad_interface_write_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: json_get_2d
    !
    ! DESCRIPTION:
    !> @brief Extract a 2D array of real(dp) values from a json file.
    !> @details
    !! This is a hopefully temporary workaround for a feature that is not implemented in the
    !! json-fortran package.
    !----------------------------------------------------------------------------------------------
    subroutine json_get_2d(jsonf, key, array)
        use json_module
        type(json_file) :: jsonf
        type(json_core) :: core
        type(json_value), pointer :: outer, inner
        character(len=*), intent(in) :: key
        real(dp), allocatable, intent(out) :: array(:,:)
        real(dp), allocatable :: work_vec(:)
        integer :: i, n1, n2
        logical :: is_matrix

        call core%initialize()
        call jsonf%get(key, outer)
        call core%matrix_info(outer, is_matrix, n_sets=n1, set_size=n2)
        if (.not. is_matrix) then
            write(stderr, *) "Error in json_get_2d."
            write(stderr, *) "  Value in "//key//" is not a matrix."
        end if
        allocate(array(n1, n2))
        do i = 1, n1
            call core%get_child(outer, i, inner)
            call core%get(inner, work_vec)
            array(i, :) = work_vec
        end do
    end subroutine json_get_2d


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: EElP
    !
    ! DESCRIPTION:
    !> @brief Calculate the electrostatic potential energy between two sets of point charges.
    !> @details
    !! Optionally, a cutoff distance can be passed to the function. If present, interactions
    !! between pairs of charges whose distance is greater than the cutoff distance will not be
    !! added to the final energy.
    !----------------------------------------------------------------------------------------------
    pure function eelp(geom1, geom2, q1, q2, cutoff) result (en)
        real(dp), intent(in) :: geom1(:, :)
        real(dp), intent(in) :: geom2(:, :)
        real(dp), intent(in) :: q1(:)
        real(dp), intent(in) :: q2(:)
        real(dp), optional, intent(in) :: cutoff
        real(dp) :: en
        integer :: i
        integer :: j
        real(dp) :: dxyz(size(geom1, 1))
        real(dp) :: d

        en = 0.0_dp
        do i = 1, size(q1)
            do j = 1, size(q2)
                dxyz = geom1(:, i) - geom2(:, j)
                d = sqrt(dot_product(dxyz, dxyz))
                if (present(cutoff)) then
                    if (d > cutoff) cycle
                end if
                en = en + q1(i) * q2(j) / d
            end do
        end do
    end function eelp

end module interface_mod
