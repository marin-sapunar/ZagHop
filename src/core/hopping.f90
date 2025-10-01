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
        use system_type_mod
        use control_var
        use decoherence_mod
        use phase_mod
        use sh_ldiab_mod
        use sh_fssh_mod
        use sh_lz_mod
        use tdc_mod
        use constants
        integer :: i
        real(dp) :: t0

        if ((ctrl%adt) .and. (tr1%step /= 0)) then
            do i = 1, tr1%wf%n_state_group
                call adt2overlap(tr2%wf%overlap(2*i)%c, tr1%wf%overlap(2*i)%c, tr1%wf%overlap(2*i-1)%c)
            end do
        end if
        !> @todo Move this to more appropriate place.
        if ((ctrl%tdc_type /= 'nadvec') .and. (ctrl%vrescale == 3)) then
            tr1%wf%need_nadv = .false.
        end if


        ctrl%hop = .false.
        select case(ctrl%sh)
        case(1)
            call lzsh(ctrl%dt, ctrl%qm_en_err, ctrl%lz_prob_conv, ctrl%lz_min_dt, &
            &         ctrl%dt_0, ctrl%rng)
        case(2)
            call decoherence()
            call phasematch()
            t0 = trajectory_data(index_offset(data_index_1, -2))%time !< @todo Clean up this type of indexing.
            call sh_adiabatic(ctrl%tdc_type, ctrl%ene_interpolate, ctrl%tdc_interpolate, t0, &
            &                 tr2%time, tr1%time, tr2%wf, tr1%wf, ctrl%shnstep, &
            &                 tr2%velo(:, tr2%qind), tr1%velo(:, tr1%qind), ctrl%rng)
        case(3)
            call decoherence()
            call phasematch()
            call sh_diabatic(tr2%time, tr1%time, tr2%wf, tr1%wf, ctrl%rng)
        end select

        if (tr2%wf%active_state /= tr1%wf%active_state) then
            ctrl%hop = .true.
            !> @todo Move this to a more appropriate place.
            ! If vrescale=3 (rescale along nonadiabatic coupling vector), ensure that the
            ! nonadiabatic coupling vector between the previous and current state is allocated
            if (ctrl%vrescale == 3) then
                tr1%wf%need_nadv(tr1%wf%active_state, tr2%wf%active_state) = .true.
                tr1%wf%need_nadv(tr2%wf%active_state, tr1%wf%active_state) = .true.
            end if
        end if
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
    subroutine sh_rescalevelo(opt_mc, opt_fh, amask, pst, wf, mass, velo)
        use system_type_mod, only : ekin
        use mqc_wave_function_mod, only : mqc_wave_function
        integer, intent(in) :: opt_mc !< Type of momentum correction.
        integer, intent(in) :: opt_fh !< Behaviour at frustrated hop.
        integer, intent(in) :: amask(:) !< Atoms considered when rescaling.
        integer, intent(in) :: pst !< Previous state.
        type(mqc_wave_function), intent(in) :: wf !< Wave function at previous step.
        real(dp), intent(in) :: mass(:) !< Masses.
        real(dp), intent(inout) :: velo(:, :) !< Velocities.
        integer :: cst !< Current state.
        real(dp) :: mvel(size(velo, 1), size(amask)) !< Mass weighted velocity.
        real(dp) :: pgrd(size(velo, 1), size(amask)) !< Gradient of the previous state.
        real(dp) :: cgrd(size(velo, 1), size(amask)) !< Gradient of the current state.
        real(dp) :: m(size(velo, 1), size(amask)) !< Temporary mass array.
        real(dp) :: rescale_dir(size(velo, 1), size(amask)) !< Direction along which to rescale.
        real(dp) :: nadv(size(velo, 1)*size(amask)) !< Nonadiabatic coupling vector between the two states.
        real(dp) :: poten(size(wf%qm_state)) !< Potential energies of all states.
        real(dp) :: mvel_dir !< Component of mass weighted velocity along rescale direction.
        real(dp) :: delta_e !< Required change in kinetic energy.

        if (stdp2) write(stdout, '(5x,a)') 'Ensuring energy conservation.'

        cst = wf%active_state
        pgrd = wf%qm_state(pst)%gradient
        cgrd = wf%qm_state(cst)%gradient
        if (opt_mc == 3) then
            nadv = wf%qm_state(pst)%nadv(cst)%c
        end if
        poten = wf%qm_state(:)%energy

        ! Work with temporary arrays and use mass-weighted coordinates.
        m = spread(mass(amask), 1, size(velo, 1))
        mvel = velo(:, amask) * sqrt(m)

        select case(opt_mc)
        case(0) ! No rescaling.
            return
        case(1) ! Rescale along velocity vector.
            rescale_dir = mvel
        case(2) ! Rescale along gradient difference vector.
            rescale_dir = (pgrd - cgrd) / sqrt(m)
        case(3) ! Rescale along nonadiabatic coupling vector.
            rescale_dir = reshape(nadv, shape(mvel)) / sqrt(m)
        end select

        ! Rescale velocity
        delta_e = poten(cst) - poten(pst)
        rescale_dir = rescale_dir / sqrt(sum(rescale_dir**2))
        mvel_dir = sum(rescale_dir * mvel)
        if (mvel_dir**2 > 2 * delta_e) then
            ! The sign here is based on the work of Herman (see chapter 6 of 10.1007/0-306-46949-9
            ! and references within). Choosing the other direction would also work, but this one
            ! ensures that the change in momentum is as small as possible (but also means that the
            ! direction of the change doesn't depend on the sign of the vector along which rescaling
            ! is performed).
            mvel = mvel + rescale_dir * (sign(sqrt(mvel_dir**2 - 2*delta_e), mvel_dir) - mvel_dir)
        else
            if (stdp1) then
                write(stdout, '(7x,a)') 'Insufficient energy for hop.'
                write(stdout, '(7x,a,e16.8)') 'Energy difference:', delta_e
                write(stdout, '(7x,a,e16.8)') 'Available energy:', mvel_dir**2 / 2
            end if

            cgrd = pgrd
            cst = pst
            select case(opt_fh)
            case(1) ! Just return to previous state.
                continue
            case(2) ! Invert velocity along rescale direction.
                mvel = mvel - 2 * mvel_dir * rescale_dir
            end select
        end if

        ! Save changes to velocity array.
        velo(:, amask) = mvel / sqrt(m)
    end subroutine sh_rescalevelo


end module hopping_mod
