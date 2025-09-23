!--------------------------------------------------------------------------------------------------
! MODULE: state_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August 2025
!
! DESCRIPTION:
!> @brief Defines type for holding all variables connected to states of the system
!--------------------------------------------------------------------------------------------------
module state_mod
    ! Import variables
    use global_defs
    implicit none

    private
    public :: state


    !----------------------------------------------------------------------------------------------
    ! TYPE: state
    !> @brief Holds all information on a single electronic state of the system.
    !----------------------------------------------------------------------------------------------
    type state
        character(len=:), allocatable :: label
        integer :: group = 0 !< Index of the group this state belongs to.
        integer :: group_state = 0 !< Index of this state within its group.
        integer :: i_state = 0 !< Index of this state in the overall state list.
        integer :: n_state = 0 !< Total number of states in the system.
        integer :: multiplicity = 1
        real(dp) :: energy
        logical :: need_gradient = .false.
        real(dp), allocatable :: gradient(:, :) !< Gradient of the state.
        type(rvec), allocatable :: nadv(:) !< Non-adiabatic vectors involving this state.
                                           !! Dimension: n_state
        real(dp), allocatable :: soc(:) !< Spin-orbit couplings involving this state.
                                        !! Dimension: n_state
    contains
        procedure :: initialize
    end type state

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: initialize
    !
    ! DESCRIPTION:
    !> @brief Allocates arrays in the state variable and initializes some variables.
    !----------------------------------------------------------------------------------------------
    subroutine initialize(self, group, group_state, i_state, n_state, multiplicity)
        class(state), intent(inout) :: self
        integer, intent(in) :: group
        integer, intent(in) :: group_state
        integer, intent(in) :: i_state
        integer, intent(in) :: n_state
        integer, intent(in) :: multiplicity

        self%group = group
        self%group_state = group_state
        self%i_state = i_state
        self%n_state = n_state
        self%multiplicity = multiplicity
        allocate(self%nadv(n_state))
    end subroutine initialize


end module state_mod