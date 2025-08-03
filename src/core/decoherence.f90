!--------------------------------------------------------------------------------------------------
! MODULE: decoherence_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
! DESCRIPTION:
!> @brief Subroutines regarding decoherence in surface hopping calculations.
!--------------------------------------------------------------------------------------------------
module decoherence_mod
    ! Import variables
    use global_defs
    implicit none

    private
    public :: decoherence

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Decoherence
    !
    ! DESCRIPTION:
    !> @brief Apply decoherence correction to electronic wave function.
    !> @details
    !! Method of decoherence correction is chosen based on the value of the ctrl%decohlvl variable.
    !----------------------------------------------------------------------------------------------
    subroutine decoherence()
        use system_var
        use control_var


        select case(ctrl%decohlvl)
        case(0) ! No decoherence correction.
        case(1) ! Energy based decoherence.
            if (stdp3) write(stdout, *) ' Dechoerence correction type EDC.'
            call edc(tr1%cstate, tr1%qkine(), tr1%qe, ctrl%dt, tr1%cwf)
        end select
    end subroutine decoherence


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: EDC
    !
    ! DESCRIPTION:
    !> @brief Energy based decoherence scheme.
    !> @details
    !! Decoherence algorithm suggested by Grannuci and Persico based on the prescription of Truhlar
    !! and co-workers in the framework of mean field methods.
    !!
    !! Constant C set to 0.1 Hartree
    !! References:
    !! - G. Granucci, M. Persico: J. Chem. Phys. 126 (2007) 134114
    !! - - DOI: 10.1063/1.2715585.
    !! - C. Zhu, A. W. Jasper, and D. G. Truhlar, J. Chem. Theory Comput. 1, 527 (2005)
    !! - - DOI: 10.1021/ct050021p
    !----------------------------------------------------------------------------------------------
    subroutine edc(cstate, kinen, stateen, dt, cwf)
        integer, intent(in)  :: cstate
        real(dp), intent(in) :: kinen
        real(dp), intent(in) :: stateen(:)
        real(dp), intent(in) :: dt
        complex(dp), intent(inout):: cwf(:)
        integer :: i
        real(dp) :: tau
        real(dp) :: sumabs
        real(dp), parameter :: Cconst = 0.1_dp
        complex(dp) :: c_temp
        real(dp), parameter :: tinydp = 1.e-15_dp

        if (kinen < tinydp) return !< To avoid division when initial velocity is zero.
        sumabs = 0.0_dp
        do i = 1, size(cwf)
            if (i /= cstate) then ! Correction for other states.
                tau = (1.0_dp + Cconst / kinen) / abs(stateen(i) - stateen(cstate))
                c_temp = cwf(i) * sqrt(exp(-1.0_dp * dt / tau))
                sumabs = sumabs + abs(c_temp)**2
                cwf(i) = c_temp
            end if
        end do
        c_temp = cwf(cstate) * sqrt(1 - sumabs) / abs(cwf(cstate))
        cwf(cstate) = c_temp ! Correction for current state.
    end subroutine edc


end module decoherence_mod
