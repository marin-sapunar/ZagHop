!--------------------------------------------------------------------------------------------------
! MODULE: tdc_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date December, 2017
!> @author Cristina Sanz, Autonoma University Madrid
!> @date May, 2024: sovec2tdc added for the time derivative couplings using spin-orbit couplings
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
    public :: adt2overlap
    public :: overlap2tdc
    public :: nadvec2tdc, sovec2tdc

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: adt2overlap
    !
    ! DESCRIPTION:
    !> @brief Calculate the overlap matrix from the adiabatic-to-diabaitc transform matrix.
    !----------------------------------------------------------------------------------------------
    subroutine adt2overlap(adt0, adt1, olap)
        use linalg_wrapper_mod
        real(dp), intent(in) :: adt0(:, :) !< A2D transform matrix at previous step.
        real(dp), intent(in) :: adt1(:, :) !< A2D transform matrix at current step.
        real(dp), allocatable, intent(out) :: olap(:, :) !< Adiabatic WF overlap matrix.

        if (.not. allocated(olap)) allocate(olap, mold=adt1)
        call gemm(adt0, adt1, olap, transa='T')
    end subroutine adt2overlap


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
    subroutine nadvec2tdc(nadvec, velo, tdc)
        real(dp), intent(in) :: nadvec(:, :, :) !< Nonadiabatic coupling vectors.
        real(dp), intent(in) :: velo(:, :) !< Velocities.
        real(dp), intent(out) :: tdc(:, :) !< Time-derivative couplings.
        real(dp), allocatable :: velo_vec(:)
        integer :: i
        integer :: j

        !allocate(tdc(size(nadvec, 2), size(nadvec, 3)))
        velo_vec = reshape(velo, [size(velo, 1)*size(velo, 2)])
        do i = 1, size(nadvec, 2)
            do j = 1, size(nadvec, 3)
                tdc(i, j) = dot_product(nadvec(:, i, j), velo_vec)
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
