module evaluator_base_mod
    use global_defs
    implicit none

    private
    public :: potential_evaluator

    type, abstract :: potential_evaluator
        integer :: n_state = 0
        integer :: n_group = 0
        integer, allocatable :: group_nstate(:)
        integer, allocatable :: group_i0(:)
    contains
        procedure(eval), deferred :: eval
        procedure(update_geometry), deferred :: update_geometry
        procedure(get_Hamiltonian), deferred :: get_Hamiltonian
        procedure(get_transformation), deferred :: get_transformation
        procedure(get_energy_single), deferred :: get_energy_single
        procedure(get_energy_all), deferred :: get_energy_all
        procedure(get_gradient), deferred :: get_gradient
        procedure(get_nadv), deferred :: get_nadv
        procedure :: set_group_nstate
        generic :: get_energy => get_energy_single, get_energy_all
    end type potential_evaluator

    
    abstract interface
        subroutine eval(self)
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
        end subroutine eval
    end interface


    abstract interface
        pure function get_energy_single(self, state) result(energy)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            integer, intent(in) :: state
            real(dp) :: energy
        end function get_energy_single
    end interface


    abstract interface
        pure function get_energy_all(self) result(energy)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            real(dp), allocatable :: energy(:)
        end function get_energy_all
    end interface


    abstract interface
        subroutine update_geometry(self, geometry)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
            real(dp), intent(in) :: geometry(:, :)
        end subroutine update_geometry
    end interface


    abstract interface
        function get_Hamiltonian(self, basis) result(H)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            real(dp), allocatable :: H(:, :)
        end function get_Hamiltonian
    end interface


    abstract interface
        function get_transformation(self, from_basis, to_basis) result(trans)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: from_basis, to_basis
            real(dp), allocatable :: trans(:, :)
        end function get_transformation
    end interface


    abstract interface
        function get_gradient(self, basis, istate) result(grad)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate
            real(dp), allocatable :: grad(:, :)
        end function get_gradient
    end interface


    abstract interface
        function get_nadv(self, basis, istate1, istate2) result(nadv)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate1, istate2
            real(dp), allocatable :: nadv(:)
        end function get_nadv
    end interface


contains


    subroutine set_group_nstate(self, nstate, multiplicity)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:)
        integer, intent(in) :: multiplicity(:)
        integer :: i, j, ci, cgroup

        self%n_group = sum(multiplicity, mask=(nstate > 0))
        allocate(self%group_nstate(self%n_group))
        allocate(self%group_i0(self%n_group))
        cgroup = 1
        ci = 1
        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            do j = 1, multiplicity(i)
                self%group_nstate(cgroup) = nstate(i)
                self%group_i0(cgroup) = ci
                cgroup = cgroup + 1
                ci = ci + nstate(i)
            end do
        end do
    end subroutine set_group_nstate


end module evaluator_base_mod