module ode_call_mod
    use global_defs
    implicit none

    private
    public :: callode

    integer :: odens !< Number of states for the ODE function.
    complex(dp), allocatable :: odecmat(:,:) !< Couplings for the ODE function.

    contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: CallODE
    !
    ! DESCRIPTION:
    !> @brief Interface for the Shampine/Gordon ODE solver.
    !----------------------------------------------------------------------------------------------
    subroutine callode(cmat, cwf, tt, edt)
        complex(dp), intent(in) :: cmat(:,:) !< Coupling matrix for the ODE function.
        complex(dp), intent(inout) :: cwf(:) !< Wave function.
        real(dp), intent(inout) :: tt !< Current time.
        real(dp), intent(in) :: edt !< Time step.

        integer :: de_flag !< ODE solver status flag.
        integer :: ns !< Number of states.
        real(dp), parameter :: relerr = 1.0d-12
        real(dp), parameter :: abserr = 1.0d-15
        real(dp), allocatable, save :: work(:)
        integer, save :: iwork(5)
        real(dp), allocatable :: yin(:)
        logical :: de_failed
        integer :: i

        odens = size(cmat, 1)
        odecmat = cmat
        if (.not. allocated(work)) allocate(work(1:100 + 42*odens))
        if (.not. allocated(yin)) allocate(yin(2 * odens))

        de_failed = .false.
        do i = 1, odens
            yin(2*i-1) = real(cwf(i))
            yin(2*i) = aimag(cwf(i))
        end do

 333    call ode(sgrhs, 2 * odens, yin, tt, tt + edt, relerr, abserr, de_flag, work, iwork)
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

 222    do i = 1, odens
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
        integer :: i, j

        y_der = 0.0_dp
        do j = 1, odens
            do i = 1, odens
                y_der(2*j - 1) = y_der(2*j - 1) + real(odecmat(j, i)) * y(2*i - 1)
                y_der(2*j - 1) = y_der(2*j - 1) - aimag(odecmat(j, i)) * y(2*i)
                y_der(2*j) = y_der(2*j) + aimag(odecmat(j, i)) * y(2*i - 1)
                y_der(2*j) = y_der(2*j) + real(odecmat(j, i)) * y(2*i)
            end do
        end do

    end subroutine sgrhs

end module ode_call_mod
