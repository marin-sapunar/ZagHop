!--------------------------------------------------------------------------------------------------
! MODULE: sh_fssh_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!> @author Cristina Sanz, Autonoma University Madrid
!> @date May, 2024: sh_sosh and sh_interpolate_sovec added for the Fewest Switches Surface 
!! Hopping method using spin-orbit couplings
!
! DESCRIPTION:
!> @brief Surface hopping algorithm.
!--------------------------------------------------------------------------------------------------
module sh_fssh_mod
    use global_defs
    implicit none

    private
    public :: sh_adiabatic, sh_sosh


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Adiabatic
    !
    ! DESCRIPTION:
    !> @brief Tully's Fewest Switches Surface Hopping Method.
    !> @details
    !! Propagates electronic wave function coefficients and determines hops for the FSSH method.
    !! The coefficients are propagated using an external subroutine for solving the system of
    !! ordinary differential equations.
    !! The nuclear time step (dt) is split into smaller time steps for which the coefficients and
    !! hopping probabilities are calculated.
    !----------------------------------------------------------------------------------------------
    subroutine sh_adiabatic(opt_clvl, opt_inte, opt_into, opt_intv, dt, nstep, qe1, qe2, cwf,   &
    &                         tst, olp1, olp2, nadv1, nadv2, vel1, vel2, fprob)
        use ode_call_mod ! Interface to Shampine/Gordon ODE solver.
        integer, intent(in) :: opt_clvl !< Method for calculating time-derivative couplings.
        integer, intent(in) :: opt_inte !< Method for interpolating energies during the time step.
        integer, intent(in) :: opt_into !< Method for interpolating overlaps during the time step.
        integer, intent(in) :: opt_intv !< Method for interpolating nad vecs during the time step.
        real(dp), intent(in) :: dt !< Nuclear dynamics time step.
        integer, intent(in) :: nstep !< Number of substeps.
        real(dp), intent(in) :: qe1(:) !< Energies at t0.
        real(dp), intent(in) :: qe2(:) !< Energies at t0 + dt.
        complex(dp), intent(inout) :: cwf(:) !< WF coefficients.
        integer, intent(inout) :: tst !< Current state.
        real(dp), intent(in) :: olp1(:, :) !< Overlap matrix between wfs at t0 - dt and t0.
        real(dp), intent(in) :: olp2(:, :) !< Overlap matrix between wfs at t0 and t0 + dt.
        real(dp), intent(in) :: nadv1(:, :, :) !< Nonadiabatic coupling vectors at t0.
        real(dp), intent(in) :: nadv2(:, :, :) !< Nonadiabatic coupling vectors at t0 + dt.
        real(dp), intent(in) :: vel1(:, :) !< Velocities at t0.
        real(dp), intent(in) :: vel2(:, :) !< Velocities at t0 + dt.
        real(dp), intent(out) :: fprob(:) !< Final probability for each state.
        real(dp) :: edt !< Time step for the propagation of the electronic WF.
        real(dp) :: tt !< Current time during propagation.
        real(dp) :: prob !< Probability of hopping into a state.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.

        integer :: i
        integer :: st
        integer :: de_flag
        real(dp) :: rnum !< Random number for surface hopping.

        odens = size(qe1)
        allocate(odeen(odens))
        allocate(odecmat(odens, odens))

        ! Propagation time step.
        tt = 0.0_dp
        edt = dt / nstep
        
        do i = 1, nstep
            fprob = 0.0_dp
            ! Get energies and TDCs for current substep.
            call sh_interpolate_energy(opt_inte, nstep, i, qe1, qe2, odeen)
            select case(opt_clvl)
            case(1)
                call sh_interpolate_overlap(opt_into, nstep, i, dt, olp1, olp2, odecmat)
            case(2)
                call sh_interpolate_nadvec(opt_intv, nstep, i, dt, nadv1, nadv2, vel1, vel2, odecmat)
            end select

            ! Propagate wf coefficients.
            call callode(odens, cwf, tt, edt, de_flag)

            ! Determine hopping probabilities.
            call random_number(rnum)
            cprob = 0.0_dp
            hop: do st = 1, odens
                if (st == tst) cycle
                prob = - 2 * edt * odecmat(st, tst) * real(conjg(cwf(st)) * cwf(tst)) / &
                     & (abs(cwf(tst))**2)
                if (prob > 0.0_dp) then ! Not actual probability, can be negative.
                    cprob = cprob + prob
                    fprob(st) = fprob(st) + prob
                    if (rnum < cprob) then
                        tst = st
                        exit hop
                    end if
                end if
            end do hop
            
        end do

        deallocate(odeen)
        deallocate(odecmat)
    end subroutine sh_adiabatic

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_sosh
    !
    ! DESCRIPTION:
    !> @brief Tully's Fewest Switches Surface Hopping Method.
    !> @details
    !! Propagates electronic wave function coefficients and determines hops for the FSSH method.
    !! The coefficients are propagated using an external subroutine for solving the system of
    !! ordinary differential equations.
    !! The nuclear time step (dt) is split into smaller time steps for which the coefficients and
    !! hopping probabilities are calculated.
    !----------------------------------------------------------------------------------------------
    subroutine sh_sosh(opt_clvl, opt_inte, opt_into, opt_intv, dt, nstep, qe1, qe2, cwf, tst, &
    &      olp1, olp2, nadv1, nadv2, vel1, vel2, sov1, sov2, spinst, fprob)
        use ode_call_mod ! Interface to Shampine/Gordon ODE solver.
        integer, intent(in) :: opt_clvl !< Method for calculating time-derivative couplings.
        integer, intent(in) :: opt_inte !< Method for interpolating energies during the time step.
        integer, intent(in) :: opt_into !< Method for interpolating overlaps during the time step.
        integer, intent(in) :: opt_intv !< Method for interpolating nad and soc vecs during the time step.
        real(dp), intent(in) :: dt !< Nuclear dynamics time step.
        integer, intent(in) :: nstep !< Number of substeps.
        real(dp), intent(in) :: qe1(:) !< Energies at t0.
        real(dp), intent(in) :: qe2(:) !< Energies at t0 + dt.
        complex(dp), intent(inout) :: cwf(:) !< WF coefficients.
        integer, intent(inout) :: tst !< Current state.
        real(dp), intent(in) :: olp1(:, :) !< Overlap matrix between wfs at t0 - dt and t0.
        real(dp), intent(in) :: olp2(:, :) !< Overlap matrix between wfs at t0 and t0 + dt.
        real(dp), intent(in) :: nadv1(:, :, :) !< Nonadiabatic coupling vectors at t0.
        real(dp), intent(in) :: nadv2(:, :, :) !< Nonadiabatic coupling vectors at t0 + dt.
        real(dp), intent(in) :: vel1(:, :) !< Velocities at t0.
        real(dp), intent(in) :: vel2(:, :) !< Velocities at t0 + dt.
        real(dp), intent(in) :: sov1(:, :) !< spin-orbit coupling vectors at t0.
        real(dp), intent(in) :: sov2(:, :) !< spin-orbit coupling vectors at t0 + dt.
        real(dp), intent(out) :: fprob(:) !< Final probability for each state.
        real(dp) :: edt !< Time step for the propagation of the electronic WF.
        real(dp) :: tt !< Current time during propagation.
        real(dp) :: prob !< Probability of hopping into a state.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.
        integer, intent(in) :: spinst(:) !< Multiplicity of each state

        integer :: i
        integer :: st,st2
        integer :: de_flag
        real(dp) :: rnum !< Random number for surface hopping.

        odens = size(qe1)
        allocate(odeen(odens))
        allocate(odecmat(odens, odens))

        ! Propagation time step.
        tt = 0.0_dp
        edt = dt / nstep
       
        do i = 1, nstep
             fprob = 0.0_dp
            ! Get energies and TDCs for current substep.
            call sh_interpolate_energy(opt_inte, nstep, i, qe1, qe2, odeen)
            
            !Checking the type of hopping, non-adiabatic or spin-orbit
            hop: do st = 1, odens
               if(st == tst) cycle 
               if ((spinst(st) == spinst(tst))) then !NAC
                  select case(opt_clvl)
                  case(1)
                     call sh_interpolate_overlap(opt_into, nstep, i, dt, olp1, olp2, odecmat)
                  case(2)
                     call sh_interpolate_nadvec(opt_intv, nstep, i, dt, nadv1, nadv2, vel1, vel2, odecmat)
                  end select
                  ! Propagate wf coefficients.
                  call callode(odens, cwf, tt, edt, de_flag)

                  ! Determine hopping probabilities.
                 call random_number(rnum)
                 cprob = 0.0_dp
                   prob = - 2 * edt * odecmat(st, tst) * real(conjg(cwf(st)) * cwf(tst)) / &
                     & (abs(cwf(tst))**2)
                 if (prob > 0.0_dp) then ! Not actual probability, can be negative.
                    cprob = cprob + prob
                    fprob(st) = fprob(st) + prob
                    if (rnum < cprob) then
                        tst = st
                        exit hop
                    end if
                 end if
                  
               else  !spin-orbit
                   call sh_interpolate_sovec(opt_intv, nstep, i, dt, sov1, sov2, odecmat)
                  
                   ! Propagate wf coefficients.
                   call callode(odens, cwf, tt, edt, de_flag)                   
                   
                   ! Determine hopping probabilities.
                   call random_number(rnum)
                   cprob = 0.0_dp
                   prob = - 2 * edt * aimag( odecmat(tst, st) * (conjg(cwf(tst)) * cwf(st)) ) / &
                   & (abs(cwf(tst))**2)
                   if (prob > 0.0_dp) then ! Not actual probability, can be negative.
                      cprob = cprob + prob
                      fprob(st) = fprob(st) + prob
                      if (rnum < cprob) then
                         tst = st
                         exit hop
                      end if
                    end if
               endif               
            end do hop
        enddo

        deallocate(odeen)
        deallocate(odecmat)
    end subroutine sh_sosh

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Interpolate_Energy
    !
    ! DESCRIPTION:
    !> @brief Energy interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns electronic state energies during the propagation of wave function coefficients in a
    !! surface hopping calculation. The interpolation is done based on the calculated energies at
    !! the beginning and end of the time nuclear time step, E(t) and E(t + dt).
    !!
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 0 - At all substeps, the energy is constant and equal to the E(t + dt).
    !! - 1 - The energies are equal to E(t) until t + dt/2, and to E(t + dt) afterwards.
    !! - 2 - The energies are a linear interpolation between E(t) and E(t + dt).
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_energy(opt, ni, i, en1, en2, ien)
        integer, intent(in) :: opt !< Type of interpolation.
        integer, intent(in) :: ni !< Number of electronic time steps.
        integer, intent(in) :: i !< Current electronic time step.
        real(dp), intent(in) :: en1(:) !< Energies in previous nuclear time step.
        real(dp), intent(in) :: en2(:) !< Energies in current nuclear time step.
        real(dp), intent(out) :: ien(:) !< Interpolated energies.

        select case(opt)
        ! Constant value.
        case(0)
            if (i == 1) ien = en2
        ! Constant with step in middle.
        case(1)
            if (i < ni / 2) then
                ien = en1
            else
                ien = en2
            end if
        ! Linear interpolation.
        case(2)
            ien = en1 + (en2 - en1) * (i - 1) / ni
        end select
    end subroutine sh_interpolate_energy


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Interpolate_Overlap
    !
    ! DESCRIPTION:
    !> @brief TDC interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns time-derivative couplings during the propagation of wave function coefficients in a
    !! surface hopping calculation. The interpolation is done based on the overlaps of the wfs at
    !! the beginning and end of the nuclear time step.
    !!
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 1 - Finite differences method to calculate tdc at \f$ t_0 + dt/2 \f$.
    !! - 2 - Finite differences method to calculate tdc at \f$ t_0 - dt/2 \f$ and \f$ t_0 + dt/2 \f$
    !!       and then use linear inter(extra)polation between \f$ t_0 \f$ and \f$ t_0 + dt \f$.
    !! - 4 - Norm-preserving interpolation of Meek and Levine. (10.1021/jz5009449).
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_overlap(opt, ni, i, dt, olp1, olp2, itdc)
        use tdc_mod
        integer, intent(in) :: opt !< Type of interpolation.
        integer, intent(in) :: ni !< Number of electronic time steps.
        integer, intent(in) :: i !< Current electronic time step.
        real(dp), intent(in) :: dt !< Nuclear time step.
        real(dp), intent(in) :: olp1(:, :) !< Overlaps between wfs at t0 - dt and t0.
        real(dp), intent(in) :: olp2(:, :) !< Overlaps between wfs at t0 and t0 + dt.
        real(dp), allocatable, intent(inout) :: itdc(:, :) !< Time-derivative couplings.
        real(dp), allocatable, save :: prev(:, :) !< TDCs at t0 - dt/2.
        real(dp), allocatable, save :: crnt(:, :) !< TDCs at t0 + dt/2.

        
        select case(opt)
        ! Constant value.
        case(1)
            if (i == 1) call overlap2tdc(dt, olp2, itdc)
        ! Linear interpolation using coupling matrix.
        case(2)
            if (i == 1) then
                call overlap2tdc(dt, olp1, prev)
                call overlap2tdc(dt, olp2, crnt)
            end if
            itdc = prev + (crnt - prev) * (0.5_dp + (i - 1) / ni)
        ! Norm-preserving interpolation (based on 10.1021/jz5009449).
        case(4)
            if (i == 1) call npi_tdc_integrated(dt, olp2, itdc)
        end select
    end subroutine sh_interpolate_overlap
    

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Interpolate_NadVec
    !
    ! DESCRIPTION:
    !> @brief TDC interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns time-derivative couplings during the propagation of wave function coefficients in a
    !! surface hopping calculation. The interpolation is done based on the nonadiabatic coupling 
    !! vectors at the beginning and end of the nuclear time step.
    !!
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 1 - Constant TDCs at t0 + dt.
    !! - 2 - Linear interpolation between TDCs at t0 and at t0 + dt.
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_nadvec(opt, ni, i, dt, nadv1, nadv2, vel1, vel2, itdc)
        use tdc_mod
        integer, intent(in) :: opt !< Type of interpolation.
        integer, intent(in) :: ni !< Number of electronic time steps.
        integer, intent(in) :: i !< Current electronic time step.
        real(dp), intent(in) :: dt !< Nuclear time step.
        real(dp), intent(in) :: nadv1(:, :, :) !< Coupling vectors at t0.
        real(dp), intent(in) :: nadv2(:, :, :) !< Coupling vectors at t0 + dt.
        real(dp), intent(in) :: vel1(:, :) !< Velocities at t0.
        real(dp), intent(in) :: vel2(:, :) !< Velocities at t0 + dt.
        real(dp), intent(inout) :: itdc(:, :) !< Time-derivative couplings.
        real(dp), allocatable, save :: prev(:, :) !< TDCs at t0.
        real(dp), allocatable, save :: crnt(:, :) !< TDCs at t0 + dt.
        integer :: nstate

        
        select case(opt)
        ! Constant value.
        case(1)
            if (i == 1) call nadvec2tdc(nadv2, vel2, itdc)
        ! Linear interpolation using coupling matrix.
        case(2)
            if (.not. allocated(prev)) then
               nstate=size(nadv1,2)
               allocate(prev(nstate,nstate))
               allocate(crnt(nstate,nstate))
            endif
            if (i == 1) then
                call nadvec2tdc(nadv1, vel1, prev)
                call nadvec2tdc(nadv2, vel2, crnt)
            end if
            itdc = prev + (crnt - prev) * (i - 1) / ni
        end select
    end subroutine sh_interpolate_nadvec

 !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Interpolate_SOVec
    !
    ! DESCRIPTION:
    !> @brief TDC interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns time-derivative couplings during the propagation of wave function coefficients in a
    !! surface hopping calculation. The interpolation is done based on the nonadiabatic coupling 
    !! vectors at the beginning and end of the nuclear time step.
    !!
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 1 - Constant TDCs at t0 + dt.
    !! - 2 - Linear interpolation between TDCs at t0 and at t0 + dt.
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_sovec(opt, ni, i, dt, sov1, sov2, itdc)
        use tdc_mod
        integer, intent(in) :: opt !< Type of interpolation.
        integer, intent(in) :: ni !< Number of electronic time steps.
        integer, intent(in) :: i !< Current electronic time step.
        real(dp), intent(in) :: dt !< Nuclear time step.
        real(dp), intent(in) :: sov1(:, :) !< Coupling vectors at t0.
        real(dp), intent(in) :: sov2(:, :) !< Coupling vectors at t0 + dt.
        real(dp), intent(inout) :: itdc(:, :) !< Time-derivative couplings.
        real(dp), allocatable, save :: prev(:, :) !< TDCs at t0.
        real(dp), allocatable, save :: crnt(:, :) !< TDCs at t0 + dt.
        integer :: nstate

        
        select case(opt)
        ! Constant value.
        case(1)
            if (i == 1) call sovec2tdc(sov2, itdc)
        ! Linear interpolation using coupling matrix.
        case(2)
            if (.not. allocated(prev)) then
               nstate=size(sov1,2)
               allocate(prev(nstate,nstate))
               allocate(crnt(nstate,nstate))
            endif
            if (i == 1) then
                call sovec2tdc(sov1, prev)
                call sovec2tdc(sov2, crnt)
            end if
            itdc = prev + (crnt - prev) * (i - 1) / ni
        end select
    end subroutine sh_interpolate_sovec

 
 
end module sh_fssh_mod
