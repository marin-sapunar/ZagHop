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
    public :: sh_adiabatic!, sh_sosh


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sh_fssh
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
    subroutine sh_adiabatic(opt_clvl, t0, t1, wf_t0, wf_t1, nstep, vel_t0, vel_t1, rng)
        use ode_call_mod ! Interface to Shampine/Gordon ODE solver.
        use mqc_wave_function_mod, only : mqc_wave_function
        use random_mod, only : rng_type
        integer, intent(in) :: opt_clvl !< Method for calculating time-derivative couplings.
        real(dp), intent(in) :: t0 !< Initial time step index.
        real(dp), intent(in) :: t1 !< Final time step index.
        type(mqc_wave_function), intent(in) :: wf_t0 !< Wave function at t0.
        type(mqc_wave_function), intent(inout) :: wf_t1 !< Wave function at t0 + dt.
        integer, intent(in) :: nstep !< Number of substeps.
        real(dp), intent(in) :: vel_t0(:, :) !< Velocities at t0.
        real(dp), intent(in) :: vel_t1(:, :) !< Velocities at t0 + dt.
        class(rng_type), allocatable, intent(inout) :: rng
       ! real(dp), intent(out) :: fprob(:) !< Final probability for each state.
        real(dp) :: dt !< Nuclear dynamics time step.
        
        real(dp) :: edt !< Time step for the propagation of the electronic WF.
        real(dp) :: tt !< Current time during propagation.
        real(dp) :: prob !< Probability of hopping into a state.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.

        integer :: cstate !< Current state index.
        integer :: i
        integer :: st
        integer :: de_flag
        real(dp) :: rnum !< Random number for surface hopping.
        real(dp), allocatable :: nadv1(:, :, :) !< Nonadiabatic coupling vectors at t0.
        real(dp), allocatable :: nadv2(:, :, :) !< Nonadiabatic coupling vectors at t0 + dt.
        real(dp), allocatable :: olp1(:, :) !< Overlap matrix between wfs at t0 - dt and t0.
        real(dp), allocatable :: olp2(:, :) !< Overlap matrix between wfs at t0 and t0 + dt.

        odens = size(wf_t1%coeff)
        allocate(odeen(odens))
        allocate(odecmat(odens, odens))

        ! Propagation time step.
        tt = 0.0_dp
        dt = t1 - t0
        edt = dt / nstep
        cstate = wf_t0%active_state

        select case(opt_clvl)
        case(1)
            olp1 = wf_t0%overlap(1)%c
            olp2 = wf_t1%overlap(1)%c
        case(2)
            nadv1 = build_nadvec_matrix(wf_t0%need_nadv, wf_t0%qm_state)
            nadv2 = build_nadvec_matrix(wf_t1%need_nadv, wf_t1%qm_state)
        end select

        wf_t1%prob = 0.0_dp

        do i = 1, nstep
            ! Get energies and TDCs for current substep.
            !> @todo Re-activate interpolation options.
            call sh_interpolate_energy(2, nstep, i, wf_t0%qm_state(:)%energy, wf_t1%qm_state(:)%energy, odeen)
            select case(opt_clvl)
            case(1)
                call sh_interpolate_overlap(2, nstep, i, dt, olp1, olp2, odecmat)
            case(2)
                call sh_interpolate_nadvec(1, nstep, i, dt, nadv1, nadv2, vel_t0, vel_t1, odecmat)
            end select
            
            ! Propagate wf coefficients.
            call callode(odens, wf_t1%coeff, tt, edt, de_flag)
            tt = tt + edt

            ! Determine hopping probabilities.
            call rng%uniform(rnum)
            cprob = 0.0_dp
            hop: do st = 1, odens
                if (st == cstate) cycle
                prob = - 2 * edt * odecmat(st, cstate) * real(conjg(wf_t1%coeff(st)) * wf_t1%coeff(cstate)) / &
                     & (abs(wf_t1%coeff(cstate))**2)
                if (prob > 0.0_dp) then ! Not actual probability, can be negative.
                    cprob = cprob + prob
                    wf_t1%prob(st) = wf_t1%prob(st) + prob
                    if (rnum < cprob) then
                        cstate = st
                        exit hop
                    end if
                end if
            end do hop
            
        end do

        wf_t1%active_state = cstate

        deallocate(odeen)
        deallocate(odecmat)
    end subroutine sh_adiabatic

    function build_nadvec_matrix(need_nadv, states) result(nadvec)
        use state_mod, only : state
        logical, intent(in) :: need_nadv(:, :)
        type(state), intent(in) :: states(:)
        real(dp), allocatable :: nadvec(:, :, :)
        integer :: i, j, n_dof

        if (.not. any(need_nadv)) then
            call errstop("sh_fssh_mod", "No nonadiabatic coupling vectors requested.", 1)
        end if
        do i = 1, size(states)
            if (.not. allocated(states(i)%nadv)) then
                call errstop("sh_fssh_mod", "Nonadiabatic coupling vectors not allocated.", 1)
            end if
            do j = 1, size(states)
                if (need_nadv(i, j)) then
                    if (.not. allocated(states(i)%nadv(j)%c)) then
                        call errstop("sh_fssh_mod", "Nonadiabatic coupling vector requested but not allocated.", 1)
                    end if
                    n_dof = size(states(i)%nadv(j)%c)
                end if
            end do
        end do

        allocate(nadvec(n_dof, size(states), size(states)), source=0.0_dp)
        do i = 1, size(states)
            do j = 1, size(states)
                if (need_nadv(i, j)) nadvec(:, i, j) = states(i)%nadv(j)%c
            end do
        end do
    end function build_nadvec_matrix


!> @todo Re-check and re-activate SOSH
!     !----------------------------------------------------------------------------------------------
!     ! SUBROUTINE: SH_sosh
!     !
!     ! DESCRIPTION:
!     !> @brief Tully's Fewest Switches Surface Hopping Method.
!     !> @details
!     !! Propagates electronic wave function coefficients and determines hops for the FSSH method.
!     !! The coefficients are propagated using an external subroutine for solving the system of
!     !! ordinary differential equations.
!     !! The nuclear time step (dt) is split into smaller time steps for which the coefficients and
!     !! hopping probabilities are calculated.
!     !----------------------------------------------------------------------------------------------
!     subroutine sh_sosh(opt_clvl, opt_inte, opt_into, opt_intv, dt, nstep, qe1, qe2, cwf, tst, &
!     &      olp1, olp2, nadv1, nadv2, vel1, vel2, sov1, sov2, spinst, fprob)
!         use ode_call_mod ! Interface to Shampine/Gordon ODE solver.
!         integer, intent(in) :: opt_clvl !< Method for calculating time-derivative couplings.
!         integer, intent(in) :: opt_inte !< Method for interpolating energies during the time step.
!         integer, intent(in) :: opt_into !< Method for interpolating overlaps during the time step.
!         integer, intent(in) :: opt_intv !< Method for interpolating nad and soc vecs during the time step.
!         real(dp), intent(in) :: dt !< Nuclear dynamics time step.
!         integer, intent(in) :: nstep !< Number of substeps.
!         real(dp), intent(in) :: qe1(:) !< Energies at t0.
!         real(dp), intent(in) :: qe2(:) !< Energies at t0 + dt.
!         complex(dp), intent(inout) :: cwf(:) !< WF coefficients.
!         integer, intent(inout) :: tst !< Current state.
!         real(dp), intent(in) :: olp1(:, :) !< Overlap matrix between wfs at t0 - dt and t0.
!         real(dp), intent(in) :: olp2(:, :) !< Overlap matrix between wfs at t0 and t0 + dt.
!         real(dp), intent(in) :: nadv1(:, :, :) !< Nonadiabatic coupling vectors at t0.
!         real(dp), intent(in) :: nadv2(:, :, :) !< Nonadiabatic coupling vectors at t0 + dt.
!         real(dp), intent(in) :: vel1(:, :) !< Velocities at t0.
!         real(dp), intent(in) :: vel2(:, :) !< Velocities at t0 + dt.
!         real(dp), intent(in) :: sov1(:, :) !< spin-orbit coupling vectors at t0.
!         real(dp), intent(in) :: sov2(:, :) !< spin-orbit coupling vectors at t0 + dt.
!         real(dp), intent(out) :: fprob(:) !< Final probability for each state.
!         real(dp) :: edt !< Time step for the propagation of the electronic WF.
!         real(dp) :: tt !< Current time during propagation.
!         real(dp) :: prob !< Probability of hopping into a state.
!         real(dp) :: cprob !< Cumulative probability of hopping into any state.
!         integer, intent(in) :: spinst(:) !< Multiplicity of each state/block (ask Graham)

!         integer :: i
!         integer :: st,st2
!         integer :: de_flag
!         real(dp) :: rnum !< Random number for surface hopping.

! !        write(69,*)nadv1(:,:,:),nadv2(:,:,:),spinst(:),qe1(:),qe2(:)
! !        call flush(69)
        
!         odens = size(qe1)
!         allocate(odeen(odens))
!         allocate(odecmat(odens, odens))

!         ! Propagation time step.
!         tt = 0.0_dp
!         edt = dt / nstep
       
!         do i = 1, nstep
!              fprob = 0.0_dp
!             ! Get energies and TDCs for current substep.
!             call sh_interpolate_energy(opt_inte, nstep, i, qe1, qe2, odeen)
            
!             !Checking the type of hopping, non-adiabatic or spin-orbit
!             hop: do st = 1, odens               
!                if(st == tst) cycle hop 
!                if ((spinst(st) == spinst(tst))) then !NAC
!                   select case(opt_clvl)
!                   case(1)
!                      call sh_interpolate_overlap(opt_into, nstep, i, dt, olp1, olp2, odecmat)
!                   case(2)
!                      call sh_interpolate_nadvec(opt_intv, nstep, i, dt, nadv1, nadv2, vel1, vel2, odecmat)
!                   end select
!                   ! Propagate wf coefficients.
!                   call callode(odens, cwf, tt, edt, de_flag)

!                   ! Determine hopping probabilities.
!                  call random_number(rnum)
!                  cprob = 0.0_dp
!                    prob = - 2 * edt * odecmat(st, tst) * real(conjg(cwf(st)) * cwf(tst)) / &
!                      & (abs(cwf(tst))**2)
!                  if (prob > 0.0_dp) then ! Not actual probability, can be negative.
!                     cprob = cprob + prob
!                     fprob(st) = fprob(st) + prob
!                     if (rnum < cprob) then
!                         tst = st
!                         exit hop
!                     end if
!                  end if                  
!                else  !spin-orbit
! !                   write(69,*) "else ",st,tst,spinst(st),spinst(tst)
! !                   call flush(69)
!                    call sh_interpolate_sovec(opt_intv, nstep, i, dt, sov1, sov2, odecmat)
       
                  
!                    ! Propagatre wf coefficients.
!                    call callode(odens, cwf, tt, edt, de_flag)    
!                                write(69,*) 'INTERPOLATE',odeen,odecmat,cwf
                   
!                    ! Determine hopping probabilities.
!                    call random_number(rnum)
!                    cprob = 0.0_dp
!                    prob = - 2 * edt * aimag( odecmat(tst, st) * (conjg(cwf(tst)) * cwf(st)) ) / &
!                    & (abs(cwf(tst))**2)
!                    if (prob > 0.0_dp) then ! Not actual probability, can be negative.
!                       cprob = cprob + prob
!                       fprob(st) = fprob(st) + prob
!                       if (rnum < cprob) then
!                          tst = st
!                          exit hop
!                       end if
!                     end if
!                endif               
!             end do hop
!         enddo

!         write(69,*)"Probability calculated with SH ", prob,tst, cwf,odecmat
!         call flush(69)
        
!         deallocate(odeen)
!         deallocate(odecmat)
!     end subroutine sh_sosh

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
