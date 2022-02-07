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
        logical :: check
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
            call sh_diabatic(t(1)%nstate, ctrl%dt, t(2)%qe, t(1)%qe, t(1)%cwf, t(1)%cstate,        &
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
        integer, intent(in) :: opt !< Type of momentum correction.
        integer, intent(in) :: amask(:) !< Atoms considered when rescaling.
        real(dp), intent(in) :: mass(:) !< Masses.
        real(dp), intent(in) :: ppe !< Previous state potential energy.
        real(dp), intent(in) :: cpe !< Current state potential energy.
        real(dp), intent(in) :: pgrd(:, :) !< Gradient on previous state.
        real(dp), intent(in) :: cgrd(:, :) !< Gradient on current state.
        real(dp), intent(inout) :: velo(:, :) !< Velocities.
        integer, intent(out) :: flag !< 0 for successful hop, 1 for rejected hop.
        real(dp) :: crit !< Available energy for rescaling.
        real(dp) :: ke !< Kinetic energy.
        real(dp) :: vscale !< Scaling factor.
        real(dp) :: v(size(velo, 1), size(amask)) !< Temporary velocity array.
        real(dp) :: m(size(amask)) !< Temporary mass array.

        ! Work with temporary arrays.
        v = velo(:, amask)
        m = mass(amask)

        flag = 1
        select case(opt)
        case(0) ! No rescaling.
        case(1) ! Rescale along velocity vector.
            ke = ekin(m, v)
            crit = ke + ppe - cpe
            if (crit < 0) return
            vscale = sqrt(crit / ke)
            v = v * vscale
            flag = 0
        case(2) ! Rescale along gradient difference vector.
            !> @todo Add rescaling along gradient difference vector.
            stop 'gdiff rescale not implemented'
        case(3) ! Rescale along nonadiabatic coupling vector.
            !> @todo Add rescaling along nonadiabatic coupling vector.
            stop 'nadvec rescale not implemented'
        end select

        ! Save changes to velocity array.
        velo(:, amask) = v
    end subroutine sh_rescalevelo


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: EKin
    !> @brief Calculate the kinetic energy for a set of particles.
    !----------------------------------------------------------------------------------------------
    pure function ekin(mass, velo) result (en)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: velo(:, :)
        real(dp) :: en
        integer :: i

        en = 0.0_dp
        do i = 1, size(mass)
            en = en + mass(i) * dot_product(velo(:, i), velo(:, i)) * 0.5_dp
        end do
    end function ekin


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
