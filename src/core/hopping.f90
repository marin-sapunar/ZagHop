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
    ! SUBROUTINE: sh_rescalevelo
    !
    ! DESCRIPTION:
    !> @brief Handle momentum correction after hops
    !> @details
    !! The behavoiur of the subroutine is determined by two options, opt_mc and opt_fh.
    !! The first, opt_mc, determines the direction along which momentum should be rescaled:
    !! - 0 - No velocity rescaling (shouldn't be selected)
    !! - 1 - Rescale along the velocity vector (shouldn't be selected)
    !! - 2 - Rescale along the component along the gradient difference vector
    !! - 3 - Rescale along the component along the nonadiabatic coupling vector.
    !!
    !! The second option, opt_fh, determines what should be done for "frustrated hops" when there 
    !! is not enough kinetic energy along the chosen direction to conserve the total energy after
    !! a hop:
    !! - 1 - Return to previous state without changing velocity.
    !! - 2 - Return to previous state, but also invert the velocity along the rescale direction.
    !!       (this option only makes sense with opt_mc=2 or 3)
    !----------------------------------------------------------------------------------------------
    subroutine sh_rescalevelo(opt_mc, opt_fh, amask, pst, cst, mass, poten, pgrd, cgrd, nadv, velo)
        use system_var, only : ekin
        integer, intent(in) :: opt_mc !< Type of momentum correction.
        integer, intent(in) :: opt_fh !< Behaviour at frustrated hop.
        integer, intent(in) :: amask(:) !< Atoms considered when rescaling.
        integer, intent(in) :: pst !< Previous state.
        integer, intent(inout) :: cst !< Current state.
        real(dp), intent(in) :: mass(:) !< Masses.
        real(dp), intent(in) :: poten(:) !< Potential energies of the states.
        real(dp), intent(in) :: pgrd(:, :) !< Gradient on previous state.
        real(dp), intent(inout) :: cgrd(:, :) !< Gradient on current state.
        real(dp), intent(in) :: nadv(:, :, :) !< Nonadiabatic coupling vectors
        real(dp), intent(inout) :: velo(:, :) !< Velocities.
        real(dp) :: v(size(velo, 1), size(amask)) !< Temporary velocity array.
        real(dp) :: m(size(amask)) !< Temporary mass array.
        real(dp) :: rescale_dir(size(velo, 1), size(amask)) !< Direction along which to rescale.
        real(dp) :: absv_dir !< Speed along rescale direction.
        real(dp) :: newv_dir !< New speed along rescale direction.
        real(dp) :: de !< Required change in kinetic energy.

        ! Work with temporary arrays.
        v = velo(:, amask)
        m = mass(amask)

        select case(opt_mc)
        case(0) ! No rescaling.
            continue
        case(1) ! Rescale along velocity vector.
            rescale_dir = v
        case(2) ! Rescale along gradient difference vector.
            rescale_dir = pgrd - cgrd
        case(3) ! Rescale along nonadiabatic coupling vector.
            rescale_dir = reshape(nadv(:, pst, cst), shape(v))
        end select

        ! Rescale velocity
        de = poten(cst) - poten(pst)
        rescale_dir = rescale_dir / sqrt(sum(rescale_dir**2))
        absv_dir = sum(v * rescale_dir)
        if (de < ekin(m, absv_dir * rescale_dir)) then
            newv_dir = sqrt(absv_dir**2 - 2 * de / sum(spread(mass, 1, 3) * rescale_dir**2))
            v = v + rescale_dir * (newv_dir - absv_dir)
        else
            write(stdout, '(7x,a)') 'Insufficient energy for hop.'
            write(stdout, '(7x,a,e16.8)') 'Energy difference:', de
            write(stdout, '(7x,a,e16.8)') 'Available energy:', ekin(m, absv_dir * rescale_dir)

            cgrd = pgrd
            cst = pst
            select case(opt_fh)
            case(1) ! Just return to previous state.
                continue
            case(2) ! Invert velocity along rescale direction.
                v = v - 2 * absv_dir * v
            end select
        end if

        ! Save changes to velocity array.
        velo(:, amask) = v
    end subroutine sh_rescalevelo


end module hopping_mod
