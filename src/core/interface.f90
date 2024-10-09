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
        use system_var, only : trajtype
        use matrix_mod, only : unit_mat
        use file_mod, only : reader
        use model_mod, only : qmodel
        type(trajtype), intent(inout) :: t
        type(reader) :: readf
        logical :: hop
        integer :: cunit, i, j, d1, d2
        logical :: check(5)
        character(len=200) :: yamfmt
        real(dp), allocatable :: rn1(:)
        real(dp), allocatable :: rn2(:, :)
        real(dp), allocatable :: rn3(:, :, :)

        if (t%ndim == 1) then
            write(yamfmt, '(a,i0,a)') "('    - [ ', e22.12, ' ]')"
        else
            write(yamfmt, '(a,i0,a)') "('    - [ ',  ", t%ndim-1, "(e22.12, ' ,'), e22.12, ' ]')"
        end if

        select case (ctrl%qlib)
        case(0)
            if ((t%step > 0) .and. (.not. hop)) then
                call system('rm -rf prevstep')
                call system('cp -r '//ctrl%qmdir//' prevstep')
            end if

            open(newunit=cunit, file='qm.yaml', action='write')
            write(cunit, '(a)') 'system:'
            write(cunit, '(2x, a, i0)') "step : ", t%step
            write(cunit, '(2x, a)') "geom : "
            do i = 1, t%qnatom
                write(cunit, yamfmt) t%geom(:, t%qind(i))
            end do
            write(cunit, '(2x,a)') "states :"
            write(cunit, '(4x,a,i0)') "- nstate : ", t%nstate
            write(cunit, *)
            write(cunit, '(a)') "request : "
            write(cunit, '(2x,a)') "energy : True"
            write(cunit, '(2x,a,i0)') "gradient : ", t%cstate
            if (ctrl%oscill) then
                write(cunit, '(2x,a)') "oscillator_strength : True"
            end if
            close(cunit)

            if (ctrl%mm) then
                open(newunit=cunit, file='mm_geom', action='write')
                do i = 1, t%mnatom
                    write(cunit, *) t%geom(:, t%mind(i))
                end do
                close(cunit)
            end if

            ! Call interface.
            call system('rm -f qm_out.yaml')
            call system(ctrl%qprog)

            ! Check if energy and gradient files were created.
            inquire(file='qm_out.yaml', exist=check(1))
            if (.not. check(1)) then
                write(stderr, *) 'Error, qm_out.yaml file not found after QM calculation.'
                stop
            end if

            call readf%open('qm_out.yaml', abort_on_eof=.false.)
            check = .false.
            do
                call readf%next()
                if (is_iostat_end(readf%iostat)) exit
                call readf%parseline(' :,[]')
                select case(readf%args(1)%s)
                case('energy')
                    t%qe = 0.0_dp
                    do i = 1, t%nstate
                        call readf%next()
                        call readf%parseline(' ')
                        read(readf%args(2)%s, *) t%qe(i)
                    end do
                    check(1) = .true.
                case('gradient')
                    t%grad = 0.0_dp
                    do i = 1, t%qnatom
                        do j = 1, t%ndim
                            call readf%next()
                            call readf%parseline(' ')
                            read(readf%args(readf%narg)%s, *) t%grad(j, t%qind(i))
                        end do
                    end do
                    check(2) = .true.
                case('oscillator_strength')
                    if (.not. ctrl%oscill) continue
                    t%qo = 0.0_dp
                    do i = 1, t%nstate - 1
                        call readf%next()
                        call readf%parseline(' ')
                        read(readf%args(2)%s, *) t%qo(i)
                    end do
                    check(3) = .true.
                case default
                    continue
                end select
            end do
            if (.not. (check(1) .and. check(2))) then
                write(stderr, *) "Energy/gradient not found in QM output."
                stop
            end if

            if ((ctrl%tdc_type == 1) .and. (t%step /= 0)) then
                call system(ctrl%oprog)
                t%olap = unit_mat(t%max_nstate)
                open(newunit=cunit, file='qm_olap', action='read')
                read(cunit, *) d1, d2
                do i = 1, min(d1, d2)
                    read(cunit, *) t%olap(i, 1:min(d1, d2))
                end do
                close(cunit)

                do i = 1, t%max_nstate
                    if (.not. ctrl%couple(i)) then
                        t%olap(i, :) = 0.0_dp
                        t%olap(:, i) = 0.0_dp
                    end if
                end do
            else if (ctrl%tdc_type == 2) then
                open(newunit=cunit, file='qm_nadvec', action='read')
                do i = 1, t%max_nstate
                    do j = 1, t%max_nstate
                        read(cunit, *) t%nadv(:, i, j)
                    end do
                end do
                close(cunit)
            end if
        case(1)
#ifdef QUANTICS
            ! Assume all atoms are QM so no temporary arrays need to be created.
            call shzagreb_run(t%step, t%geom, t%cstate, t%qe, t%grad, t%nadv, t%adt)
#else
            write(stderr, *) 'Error in run_qm. Code not compiled with quantics interface.'
            stop
#endif
        case(2)
            call qmodel%eval(t%geom, t%cstate, t%qe, t%grad, t%nadv)
        case default
            write(stderr, *) 'Error in run_qm. Unrecognized QM interface.'
            stop
        end select

        ! Add random noise to evaluated values if requested.
        if (ctrl%noise > 0.0_dp) then
            allocate(rn1(t%nstate))
            allocate(rn2(t%ndim, t%natom))
            call random_number(rn1)
            call random_number(rn2)
            rn1 = (rn1 - 0.5_dp) * ctrl%noise
            rn2 = (rn2 - 0.5_dp) * ctrl%noise
            t%qe = t%qe + rn1
            t%grad = t%grad + rn2
            if (allocated(t%nadv)) then
                allocate(rn3(t%natom * t%ndim, t%nstate, t%nstate))
                call random_number(rn3)
                rn3 = (rn3 - 0.5_dp) * ctrl%noise
                t%nadv = t%nadv + rn3
            end if
        end if
    end subroutine run_qm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Interface_MMRun
    !----------------------------------------------------------------------------------------------
    subroutine run_mm(t)
        use system_var, only : trajtype
        type(trajtype), intent(inout) :: t
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
