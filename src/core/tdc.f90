!--------------------------------------------------------------------------------------------------
! MODULE: tdc_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date December, 2017
!
! DESCRIPTION: 
!> @brief Contains subroutines for calculating the time-derivative couplings.
!--------------------------------------------------------------------------------------------------
module tdc_mod
    use global_defs
    use npi_mod, only : npi_tdc_integrated
    implicit none

    private
    public :: npi_tdc_integrated
    public :: overlap2tdc
    public :: nadvec2tdc, sovec2tdc

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Overlap2TDC
    !
    ! DESCRIPTION:
    !> @brief Calculate the time-derivative couplings using the overlap matrix.
    !> @details
    !! The coupling matrix at time \f$ t_0 + dt/2 \f$
    !! \f[
    !! D_{ij} = \braket{\psi_i(t_0 + dt / 2)}{\dv{\psi_j(t_0 + dt / 2)}{t}}
    !! \f]
    !! is calculated by the finite differences method where
    !! \f[
    !! \bra{\psi_i(t_0 + dt/2)} \approx \frac{1}{2} (\bra{\psi_i(t_0 + dt)} + 
    !!                                        \bra{\psi_i(t_0})
    !! \f]
    !! and:
    !! \f[
    !! \ket{\dv{\psi_i(t_0 + dt/2)}{t}} \approx \frac{1}{2} (\bra{\psi_i(t_0 + dt)} - 
    !!                                                \bra{\psi_i(t_0})
    !! \f]
    !! so the final equation is:
    !! \f[
    !! D_{ij} \approx \frac{1}{2 dt} (\braket{\psi_i(t_0}{\psi_j(t_0 + dt)} -
    !!                                       \braket{\psi_i(t_0 + dt)}{\psi_j(t_0)}
    !! \f]
    !----------------------------------------------------------------------------------------------
    subroutine overlap2tdc(dt, olp, tdc)
        real(dp), intent(in) :: dt !< Time step.
        real(dp), intent(in) :: olp(:, :) !< Overlaps between wfs at t0 and t0+dt
        real(dp), allocatable, intent(out) :: tdc(:, :) !< Time-derivative couplings.

        if (.not. allocated(tdc)) allocate(tdc, mold=olp)
        tdc = (olp - transpose(olp)) * 0.5_dp / dt
    end subroutine overlap2tdc


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: NadVec2TDC
    !
    ! DESCRIPTION:
    !> @brief Calculate the time-derivative couplings using the nonadiabatic coupling vectors.
    !> @details
    !! The coupling between states i and j is calculated as a scalar product between the
    !! nonadiabatic coupling vector and the velocity vector.
    !----------------------------------------------------------------------------------------------
    subroutine nadvec2tdc(nadvec, velo, cmat)
        real(dp), intent(in) :: nadvec(:, :, :) !< Nonadiabatic coupling vectors.
        real(dp), intent(in) :: velo(:, :) !< Velocities.
        real(dp), intent(out) :: cmat(:, :) !< Time-derivative couplings.
        integer :: nstate
        integer :: natom
        integer :: ndim
        integer :: i
        integer :: j
        integer :: k
        integer :: d
        integer :: c

        ndim = size(velo, 1)
        natom = size(velo, 2)
        nstate = size(cmat,1)

        cmat = 0.0_dp
        do i = 1, nstate
        do j = 1, nstate
            c = 0
            do k = 1, natom
                do d = 1, ndim
                    c = c + 1
                    cmat(i, j) = cmat(i, j) + nadvec(c, i, j) * velo(d, k)
                end do
            end do
        end do
        end do
    end subroutine nadvec2tdc

!----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SOVec2TDC
    !
    ! DESCRIPTION:
    !> @brief Calculate the time-derivative couplings using the spin-orbit coupling vectors.
    !> @details
    !! The coupling between states i and j is calculated as a scalar product between the
    !! spin-orbit coupling vector.
    !----------------------------------------------------------------------------------------------
    subroutine sovec2tdc(sovec, cmat)
        real(dp), intent(in) :: sovec(:, :) !< spin-orbit coupling vectors.
        real(dp), intent(out) :: cmat(:, :) !< Time-derivative couplings.
        integer :: nstate
        integer :: i
        integer :: j
        integer :: k
        integer :: d
        integer :: c

        nstate = size(cmat,1)

        cmat = 0.0_dp
        do i = 1, nstate
        do j = 1, nstate           
           cmat(i, j) = cmat(i, j) + sovec(i, j) 
        end do
        end do
    end subroutine sovec2tdc

end module tdc_mod
