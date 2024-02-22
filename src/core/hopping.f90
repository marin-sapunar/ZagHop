!--------------------------------------------------------------------------------------------------
! MODULE: hopping_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
! DESCRIPTION:
!> @brief Surface hopping algorithm.
!--------------------------------------------------------------------------------------------------
module hopping_mod
    use global_defs
    implicit none

    private
    public :: hopping
    public :: sh_rescalevelo
    public :: sh_frustratedhop

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Hopping
    !
    ! DESCRIPTION:
    !> @brief Selection of surface hopping algorithm.
    !> @details
    !! Calls appropriate subroutine for surface hopping depending on the ctrl%sh variable. Also
    !! calls the decoherence correction subroutine and phase matching subroutine.
    !----------------------------------------------------------------------------------------------
    subroutine hopping()
        use system_var
        use control_var
        use decoherence_mod
        use phase_mod
        use sh_ldiab_mod
        use sh_fssh_mod
        use sh_lz_mod
        use constants
        integer :: tst

        tst = t(1)%cstate
        ctrl%hop = .false.
        select case(ctrl%sh)
        case(1)
            call lzsh()
        case(2)
            call decoherence()
            call phasematch()
            call sh_adiabatic(ctrl%tdc_type, ctrl%ene_interpolate, ctrl%tdc_interpolate,           &
            &                 ctrl%tdc_interpolate, ctrl%dt, ctrl%shnstep, t(2)%qe, t(1)%qe,       &
            &                 t(1)%cwf, t(1)%cstate, t(2)%olap, t(1)%olap, t(2)%nadv, t(1)%nadv,   &
            &                 t(2)%velo(:, t(2)%qind), t(1)%velo(:, t(1)%qind), t(1)%prob)
        case(3)
            call decoherence()
            call phasematch()
            call sh_diabatic(t(1)%max_nstate, ctrl%dt, t(2)%qe, t(1)%qe, t(1)%cwf, t(1)%cstate,    &
            &                t(1)%olap, t(1)%prob)
        end select

        if (tst /= t(1)%cstate) ctrl%hop = .true.
    end subroutine hopping


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_RescaleVelo
    !
    ! DESCRIPTION:
    !> @brief
    !> @details
    !----------------------------------------------------------------------------------------------
    subroutine sh_rescalevelo(opt, amask, mass, ppe, cpe, pgrd, cgrd, velo, flag)
        use system_var, only : ekin
        integer, intent(in) :: opt !< Type of momentum correction.
        integer, intent(in) :: amask(:) !< Atoms considered when rescaling.
        real(dp), intent(in) :: mass(:) !< Masses.
        real(dp), intent(in) :: ppe !< Previous state potential energy.
        real(dp), intent(in) :: cpe !< Current state potential energy.
        real(dp), intent(in) :: pgrd(:, :) !< Gradient on previous state.
        real(dp), intent(in) :: cgrd(:, :) !< Gradient on current state.
        real(dp), intent(inout) :: velo(:, :) !< Velocities.
        integer, intent(out) :: flag !< 0 for successful hop, 1 for rejected hop.
        real(dp) :: v(size(velo, 1), size(amask)) !< Temporary velocity array.
        real(dp) :: m(size(amask)) !< Temporary mass array.
        real(dp) :: rescale_dir(size(velo, 1), size(amask)) !< Direction along which to rescale.
        real(dp) :: absv_dir !< Speed along rescale direction.
        real(dp) :: newv_dir !< New speed along rescale direction.
        real(dp) :: de !< Required change in kinetic energy.

        ! Work with temporary arrays.
        v = velo(:, amask)
        m = mass(amask)

        flag = 1
        select case(opt)
        case(0) ! No rescaling.
        case(1) ! Rescale along velocity vector.
            rescale_dir = v
        case(2) ! Rescale along gradient difference vector.
            rescale_dir = pgrd - cgrd
        case(3) ! Rescale along nonadiabatic coupling vector.
            !> @todo Add rescaling along nonadiabatic coupling vector.
            stop 'nadvec rescale not implemented'
        end select

        ! Rescale velocity
        de = cpe - ppe
        rescale_dir = rescale_dir / sqrt(sum(rescale_dir**2))
        absv_dir = sum(v * rescale_dir)
        if (de < ekin(m, absv_dir * rescale_dir)) then
            flag = 0
            newv_dir = sqrt(absv_dir**2 - 2 * de / sum(spread(mass, 1, 3) * rescale_dir**2))
            v = v + rescale_dir * (newv_dir - absv_dir)

            ! Save changes to velocity array.
            velo(:, amask) = v
        else
            write(stdout, '(7x,a)') 'Insufficient energy for hop.'
            write(stdout, '(7x,a,e16.8)') 'Energy difference:', de
            write(stdout, '(7x,a,e16.8)') 'Available energy:', ekin(m, absv_dir * rescale_dir)
            return
        end if
    end subroutine sh_rescalevelo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_FrustratedHop
    !
    ! DESCRIPTION:
    !> @brief Return to previous state if a hop is rejected.
    !> @details
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 1 - Return to previous state without changing velocity.
    !! - 2 - Return to previous state, but also invert the velocity along the relevant direction.
    !----------------------------------------------------------------------------------------------
    subroutine sh_frustratedhop(opt, pgrd, cgrd, pst, cst)
        integer, intent(in) :: opt !< Method of treating frustrated hops.
        real(dp), intent(in) :: pgrd(:, :) !< Gradient on previous state. (New current gradient.)
        real(dp), intent(inout) :: cgrd(:, :) !< Gradient on current state.
        integer, intent(in) :: pst !< Previous state. (New current state.)
        integer, intent(inout) :: cst !< Current state.

        select case(opt)
        case(1) ! Just return to previous state.
            cgrd = pgrd
            cst = pst
        case(2) ! Invert velocity along rescale direction.
            !> @todo Implement velocity inversion on frustrated hop.
            stop 'Inversion on frustrated hop not implemented.'
        end select
    end subroutine sh_frustratedhop


end module hopping_mod
