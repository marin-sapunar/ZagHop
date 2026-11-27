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
    public :: run_mm


contains


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
        read(cunit, *) t%men
        close(cunit)

        open(newunit=cunit, file='mm_grad', action='read')
        do i = 1, t%natom
            read(cunit, *) t%grad(:, i)
        end do
        close(cunit)

        open(newunit=cunit, file='mm_pbc', action='read')
        read(cunit, *) t%pbcbox
        close(cunit)

        t%men(1) = t%men(1) - eelp(t%geom(:, t%qind), t%geom(:, t%mind), t%chrg(t%qind),             &
        &         t%chrg(t%mind), ctrl%mmcut)
    end subroutine run_mm


    subroutine nad_interface_write_atoms(nat, qnat, mnat, qind, mind, sym, geo, vel, pbc, box)
        integer, intent(in) :: nat, qnat, mnat
        integer, intent(in) :: qind(:), mind(:)
        character(len=*), intent(in) :: sym(:)
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
