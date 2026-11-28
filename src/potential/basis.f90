!--------------------------------------------------------------------------------------------------
! MODULE: basis_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2026
!
! DESCRIPTION:
!> @brief Names of the bases of the electronic states.
!> @details
!! The states of a potential are given in five bases, in two spaces:
!!
!!     space                                       diabatic       MCH       adiabatic
!!     full (spin components, H_sf + H_soc)        'diabatic'     'mch'     'adiabatic'
!!     spin-free (one state per multiplet, H_sf)   'diabatic_sf'  'mch_sf'  -
!!
!! - diabatic - The states in which a model gives its Hamiltonian, or the states closest to a set
!!   of reference states.
!! - mch - Eigenstates of the spin-free Hamiltonian (molecular Coulomb Hamiltonian) within each
!!   multiplicity. With their spin components, this is the MCH representation of Mai, Marquetand
!!   and Gonzalez, WIREs Comput. Mol. Sci. 8, e1370 (2018), doi: 10.1002/wcms.1370, Sec. 3.2.
!! - adiabatic - Eigenstates of the Hamiltonian with spin-orbit couplings (the diagonal
!!   representation of the reference).
!!
!! The full space has one state for each spin component of each state. The spin-free space has
!! one state for each multiplet. With effective spin states (see state_layout_mod), each multiplet
!! is a single state and the two spaces have the same states.
!!
!! Transformations exist only between bases of the same space. A quantity of the full space is
!! brought to the spin-free space by summing over the spin components (sum_components of
!! state_layout), which is not a transformation.
!--------------------------------------------------------------------------------------------------
module basis_mod
    implicit none

    private
    public :: basis_known
    public :: basis_is_spin_free
    public :: basis_full
    public :: basis_like


contains


    !> The name is one of the five bases.
    pure function basis_known(basis) result(known)
        character(len=*), intent(in) :: basis
        logical :: known

        select case(basis)
        case('diabatic', 'mch', 'adiabatic', 'diabatic_sf', 'mch_sf')
            known = .true.
        case default
            known = .false.
        end select
    end function basis_known


    !> The basis is a basis of the spin-free space ('diabatic_sf' or 'mch_sf').
    pure function basis_is_spin_free(basis) result(spin_free)
        character(len=*), intent(in) :: basis
        logical :: spin_free

        spin_free = (basis == 'diabatic_sf' .or. basis == 'mch_sf')
    end function basis_is_spin_free


    !> Basis of the full space with the same kind of states: 'mch' for 'mch_sf', 'diabatic' for
    !! 'diabatic_sf'. A basis of the full space is returned unchanged.
    pure function basis_full(basis) result(full)
        character(len=*), intent(in) :: basis
        character(len=:), allocatable :: full

        select case(basis)
        case('diabatic_sf')
            full = 'diabatic'
        case('mch_sf')
            full = 'mch'
        case default
            full = basis
        end select
    end function basis_full


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: basis_like
    !> @brief The basis with the states of basis in the space of other.
    !> @details
    !! For example, basis_like('diabatic', 'mch_sf') is 'diabatic_sf', the diabatic basis from
    !! which the 'mch_sf' basis can be reached. The spin-free space has no adiabatic basis.
    !----------------------------------------------------------------------------------------------
    pure function basis_like(basis, other) result(res)
        character(len=*), intent(in) :: basis
        character(len=*), intent(in) :: other !< Basis whose space is used.
        character(len=:), allocatable :: res

        res = basis_full(basis)
        if (.not. basis_is_spin_free(other)) return
        select case(res)
        case('diabatic')
            res = 'diabatic_sf'
        case('mch')
            res = 'mch_sf'
        case default
            error stop 'Error in basis_mod, basis_like function. No such spin-free basis.'
        end select
    end function basis_like


end module basis_mod
