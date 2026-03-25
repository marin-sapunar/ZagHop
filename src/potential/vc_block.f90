module vc_block_mod
    use global_defs
    implicit none

    private
    public :: vc_block

    
    type vc_block
        integer :: max_order = -1
        integer :: nmode = 0
        integer :: nstate = 0
        real(dp), allocatable :: zero(:, :)
        real(dp), allocatable :: linear(:, :, :)
        real(dp), allocatable :: quad(:, :, :, :)
    contains
        procedure :: allocate => vc_block_allocate
        procedure :: check_order => vc_block_check_order
        procedure :: add_zero => vc_block_zero
        procedure :: add_linear => vc_block_linear
        procedure :: eval => vc_block_eval
    end type vc_block


contains


    subroutine vc_block_allocate(self, max_order, nmode, nstate)
        class(vc_block) :: self
        integer, intent(in) :: max_order
        integer, intent(in) :: nmode
        integer, intent(in) :: nstate

        self%max_order = max_order
        self%nmode = nmode
        self%nstate = nstate
        if (max_order >= 0) allocate(self%zero(nstate, nstate), source=0.0_dp)
        if (max_order >= 1) allocate(self%linear(nmode, nstate, nstate), source=0.0_dp)
        if (max_order >= 2) allocate(self%quad(nmode, nmode, nstate, nstate), source=0.0_dp)
    end subroutine vc_block_allocate


    subroutine vc_block_check_order(self, order, nmode, nstate)
        class(vc_block) :: self
        integer, intent(in) :: order
        integer, intent(in) :: nmode
        integer, intent(in) :: nstate

        if (self%max_order < 0) then
            call self%allocate(order, nmode, nstate)
            return
        end if
        
        if (self%nmode /= nmode .or. self%nstate /= nstate) then
            write(stderr, *) 'Error in vibronic_mod, vc_block_check_order subroutine.'
            write(stderr, *) '  Block already allocated with different dimensions.'
            stop
        end if

        if (self%max_order >= order) return

        self%max_order = order
        if (order >= 1) then
            allocate(self%linear(self%nmode, self%nstate, self%nstate), source=0.0_dp)
        end if
        if (order >= 2) then
            allocate(self%quad(self%nmode, self%nmode, self%nstate, self%nstate), source=0.0_dp)
        end if   
    end subroutine vc_block_check_order


    subroutine vc_block_zero(self, nmode, zero)
        class(vc_block) :: self
        integer, intent(in) :: nmode
        real(dp), intent(in) :: zero(:, :)

        call self%check_order(0, nmode, size(zero, 2))
        self%zero = zero
    end subroutine vc_block_zero
    

    subroutine vc_block_linear(self, linear)
        class(vc_block) :: self
        real(dp), intent(in) :: linear(:, :, :)

        call self%check_order(1, size(linear, 1), size(linear, 3))
        self%linear = linear
    end subroutine vc_block_linear


    function vc_block_eval(self, q) result(val)
        class(vc_block) :: self
        real(dp), intent(in) :: q(:)
        integer :: i, j, imode, jmode
        real(dp), allocatable :: val(:, :)

        if (self%max_order < 0) then
            write(stderr, *) 'Error in vibronic_mod, vc_block_eval function.'
            write(stderr, *) '  Block not allocated.'
            stop
        end if

        if (self%max_order == 0) then
            val = self%zero
            return
        end if
        
        allocate(val(self%nstate, self%nstate), source=self%zero)
        if (self%max_order >= 1) then
            do imode = 1, self%nmode
                do i = 1, self%nstate
                    do j = 1, self%nstate
                        val(i, j) = val(i, j) + self%linear(imode, i, j) * q(imode)
                    end do
                end do
            end do
        end if
        if (self%max_order >= 2) then
            do imode = 1, self%nmode
                do jmode = 1, self%nmode
                    do i = 1, self%nstate
                        do j = 1, self%nstate
                            val(i, j) = val(i, j) + self%quad(imode, jmode, i, j) * q(imode) * q(jmode)
                        end do
                    end do
                end do
            end do
        end if
    end function vc_block_eval


end module vc_block_mod