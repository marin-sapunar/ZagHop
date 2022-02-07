module ode_call_mod
    use global_defs
    implicit none


    integer :: odens !< Number of states for the ODE function.
    real(dp), allocatable :: odeen(:) !< Energies for the ODE function.
    real(dp), allocatable :: odecmat(:,:) !< Couplings for the ODE function.

    contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: CallODE
    !
    ! DESCRIPTION:
    !> @brief Interface for the Shampine/Gordon ODE solver.
    !----------------------------------------------------------------------------------------------
    subroutine callode(ns, cwf, tt, edt, de_flag)
        integer, intent(in) :: ns !< Number of states.
        complex(dp), intent(inout) :: cwf(:) !< Wave function.
        real(dp), intent(inout) :: tt !< Current time.
        real(dp), intent(in) :: edt !< Time step.
        integer, intent(inout) :: de_flag !< ODE solver status flag.

        real(dp), parameter :: relerr = 1.0d-12
        real(dp), parameter :: abserr = 1.0d-15
        real(dp), allocatable, save :: work(:)
        integer, save :: iwork(5)
        real(dp) :: yin(2 * ns)
        logical :: de_failed
        integer :: i

        if (.not. allocated(work)) allocate(work(1:100 + 42*ns))

        de_failed = .false.
        do i = 1, ns
            yin(2*i-1) = real(cwf(i))
            yin(2*i) = aimag(cwf(i))
        end do

 333    call ode(sgrhs, 2 * ns, yin, tt, tt + edt, relerr, abserr, de_flag, work, iwork)
        select case(de_flag)
        case(2)
            goto 222
        case(3:5)
            goto 333
        case(6)
            de_flag = 1
            if (de_failed) then
                write(stderr, '(a)') 'Error in ODEMod, CallODE subroutine.'
                write(stderr, '(a)') 'Invalid input_mod to ODE subroutine.'
                stop
            else
                de_failed = .true.
                goto 333
            end if
        case default
            write(stderr, '(a)') 'Error in ODEMod, CallODE subroutine.'
            write(stderr, '(3x,a,i0)') 'Error in ODE subroutine. Flag value: ', de_flag
        end select

 222    do i = 1, ns
            cwf(i) = cmplx(yin(2*i - 1), yin(2*i), kind = dp)
        end do
        
    end subroutine callode


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SGRHS
    !
    ! DESCRIPTION:
    !> @brief Supplies the right hand side of the ODE for the Shampine/Gordon solver.
    !----------------------------------------------------------------------------------------------
    subroutine sgrhs(ttt, y, y_der)
        real(dp) :: ttt !< Time.
        real(dp), intent(in) :: y(2*odens) !< Dependent variable.
        real(dp), intent(inout) :: y_der(2*odens) !< Value of the derivative.
        integer :: i1
        integer :: i2
        complex(dp) :: ccc(odens)
        complex(dp) :: ccc_der(odens)

        ccc = cmplx((0.0_dp, 0.0_dp), kind = dp)
        ccc_der = cmplx((0.0_dp, 0.0_dp), kind = dp)

        do i1 = 1, odens
            ccc(i1) = cmplx(y(2*i1-1), y(2*i1), kind = dp)
        end do

        do i1 = 1, odens
            do i2 = 1, odens
                ccc_der(i1) = ccc_der(i1) - cmplx(odecmat(i1, i2), 0.0_dp, kind = dp) * ccc(i2)
            end do
            ccc_der(i1) = ccc_der(i1) - cmplx(0.0_dp, odeen(i1), kind = dp) * ccc(i1)
        end do

        y_der = 0.0_dp
        do i1 = 1, odens
            y_der(2*i1 - 1) = real(ccc_der(i1))
            y_der(2*i1) = aimag(ccc_der(i1))
        end do
    end subroutine sgrhs

end module ode_call_mod
