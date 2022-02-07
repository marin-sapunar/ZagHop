!--------------------------------------------------------------------------------------------------
! MODULE: NPIMod
!
!> @brief Subroutines for the norm-preserving interpolation (NPI) method of Meek and Levine.
!> @details
!! NPI is an overlap based method for interpolating wave functions and time-derivative couplings 
!! (TDCs) during a nuclear time step in nonadiabatic dynamics calculations.
!! The method is based on the interpolation of the wave functions from those at the beginning of 
!! the time step to those at the end using a time-dependent transformation matrix:
!! \f[
!! \ket{\psi_j(t)} = \vb{U}(t) \ket{\psi_j(t_0)}
!! \f]
!! where the elements of \f$ \vb{U}(t) \f$ are:
!! \f[
!! U_{jj}(t) = \cos(\frac{\cos^{-1}W_{jj}}{\Delta t} (t - t_0))
!! \f]
!! and
!! \f[
!! U_{jk}(t) = \sin(\frac{\sin^{-1}W_{jk}}{\Delta t} (t - t_0)).
!! \f]
!! The resulting TDCs are:
!! \f[
!! \braket{\psi_k(t)}{\pdv{\psi_j(t)}{t}} = 
!! \mel{\psi_k(t_0)}{\vb{U}^\dagger(t)\pdv{t}\vb{U(t)}}{\psi_j(t_0)}
!! \f]
!! 
!! These TDCs can be analytically integrated over a given time step \f$ t_0 \f$ to \f$ t_0 + 
!! \Delta t \f$ (\ref npi_tdc_integrated subroutine). This is done in the original article and the
!! value is used as a constant during the time step.
!!
!! Original article: <a href="http://www.dx.doi.org/10.1021/jz5009449">DOI: 10.1021/jz5009449</a>
! 
!> @note The main input_mod for the NPI method is the overlap matrix containing the overlaps between
!! the electronic state wave functions at \f$ t_0 \f$ and at \f$ t_0 + \Delta t \f$. It is assumed
!! that the phases of the two sets of wave functions are matched, meaning that the change from
!! \f$ \qty{\ket{\psi_i(t_0)}} \f$ to \f$ \qty{\ket{\psi_i(t_0 + \Delta t)}} \f$ is (approximately) 
!! a rotation.
!--------------------------------------------------------------------------------------------------
module npi_mod
    use global_defs
    implicit none

    private
    public :: npi_tdc_integrated


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: NPI_TDC_Integrated
    !
    ! DESCRIPTION:
    !> @brief Average value of the TDC over a time interval \p dt.
    !> @details
    !! Returns the analytic solution of the integral:
    !! \f[
    !! \eval{\braket{\psi_k(t)}{\pdv{\psi_j(t)}{t}}}_{t_0+\Delta t/2} = 
    !! \frac{\int_{t_0}^{t_0+\Delta t}
    !!       {\mel{\psi_k(t_0)}{\vb{U}^\dagger(t)\pdv{t}\vb{U(t)}}{\psi_j(t_0)} dt}}{\Delta t}
    !! \f]
    !! 
    !! Equations are taken from the Supporting Information of the 
    !! <a href="http://www.dx.doi.org/10.1021/jz5009449">article</a>.
    !! Variable names are chosen to match those in the article. The compile time parameter \p thr
    !! can be used to change the value at which the equations are replaced by their limits as the
    !! denominators approach zero.
    !----------------------------------------------------------------------------------------------
    subroutine npi_tdc_integrated(dt, ww, tdc)
        real(dp), intent(in) :: dt !< Time step.
        real(dp), intent(in) :: ww(:, :) !< Wave function overlap matrix.
        real(dp), intent(out) :: tdc(:, :) !< Approximate time derivative coupling at t0+dt/2.
        real(dp), parameter :: thr = 0.00000001_dp !< Denominator threshod below which values are
                                                   !! replaced by their limits.
        real(dp) :: w(size(ww, 1), size(ww, 2)) !< Temporary w array.
        real(dp) :: acs(size(ww, 1), size(ww, 2)) !< ArcCos of w array elements.
        real(dp) :: ass(size(ww, 1), size(ww, 2)) !< ArcSin of w array elements.
        real(dp) :: a, b, c, d, e !< Temporary values.
        real(dp) :: denom, wlj, wlk, aslj, aslk !< Temporary values.
        integer :: k, j, n

        w = ww
        where (abs(w)>1.0_dp) w=int(w)
        acs = acos(w)
        ass = asin(w)
        n = size(w, 1)

        do k = 1, n
            do j = 1, n
                if (j == k) then
                    tdc(k, j) = 0.0_dp
                    cycle
                end if
                a = -1.0_dp !< Limit for a.
                b = 1.0_dp !< Limit for b.
                c = 1.0_dp !< Limit for c.
                d = 1.0_dp !< Limit for d.
                denom = acs(j, j) - ass(j, k)
                if (abs(denom) > thr) a = -sin(denom) / denom
                denom = acs(j, j) + ass(j, k)
                if (abs(denom) > thr) b = sin(denom) / denom
                denom = acs(k, k) - ass(k, j)
                if (abs(denom) > thr) c = sin(denom) / denom
                denom = acs(k, k) + ass(k, j)
                if (abs(denom) > thr) d = sin(denom) / denom

                e = 0.0_dp !< Limit for e when wlj = 0.
                if (n > 2) then
                    wlj = sqrt(1 - w(j, j)**2 - w(k, j)**2)
                    if (wlj > thr) then
                        wlk = - (w(j, k) * w(j, j) + w(k, k) * w(k, j)) / wlj
                        aslj = asin(wlj)
                        aslk = asin(wlk)
                        denom = aslj**2 - aslk**2
                        if (abs(denom) > thr) then
                            e = 2 * aslj * ( wlj * wlk * aslj +                                    &
                            &   (sqrt((1 - wlj**2) * (1 - wlk**2)) - 1) * aslk) / denom
                        else
                            e = wlj**2 !< Limit for e when denom = 0.
                        end if
                    end if
                end if
                ! Integrated expression for the NPI TDC.
                tdc(k, j) = (acs(j, j) * (a + b) + ass(k, j) * (c + d) + e) * 0.5_dp / dt
            end do
        end do
    end subroutine npi_tdc_integrated


end module npi_mod
