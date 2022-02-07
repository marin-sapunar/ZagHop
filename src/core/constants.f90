!--------------------------------------------------------------------------------------------------
! MODULE: constants
!
! DESCRIPTION:
!> @brief Definitions of constants and unit conversion factors.
!--------------------------------------------------------------------------------------------------
module constants
    use global_defs, only : dp

    !----------------------------------------------------------------------------------------------
    ! Unit conversions.
    !----------------------------------------------------------------------------------------------
    real(dp), parameter :: Da_me = 1822.88849_dp !< Dalton (unified atomic mass) to electron rest mass.
    real(dp), parameter :: a0_A = 0.52917721067_dp !< Bohr radius to Angstrom.
    real(dp), parameter :: aut_fs = 0.02418884326509_dp !< Atomic unit of time to femtosecond.
    real(dp), parameter :: aut_ps = 0.00002418884326509_dp !< Atomic unit of time to picosecond.
    real(dp), parameter :: eh_k = 315775.13_dp !< Hartree to Kelvin.
    real(dp), parameter :: eh_cm1 = 219474.6313702_dp !< Hartree to reciprocal centimeter.
    real(dp), parameter :: eh_eV = 27.21138602_dp !< Hartree to electronvolt.
    real(dp), parameter :: kcalmol_eh = 0.001593601_dp !< kcal/mol to Hartree.

 
end module constants
