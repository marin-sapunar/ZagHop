module vc_evaluator_mod
    use potential_evaluator_mod
    use vc_model_mod
    use rmat_mod, only : rmat
    implicit none

    type, extends(potential_evaluator) :: vc_evaluator
        type(vc_model), pointer :: model => null()
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: w_full(:, :)
        real(dp), allocatable :: w_eigvec(:, :)
        real(dp), allocatable :: w_grad(:, :, :)
        real(dp), allocatable :: dm(:, :, :)
    contains
        procedure :: initialize => vc_initialize
        procedure :: update_geometry => vc_update_geometry
        procedure :: eval => vc_eval
        procedure :: get_nadv => vc_get_nadv
        procedure :: get_oscill => vc_get_oscill
    end type vc_evaluator

contains

    subroutine vc_initialize(self, model)
        class(vc_evaluator) :: self
        type(vc_model), intent(in), target :: model

        self%model => model
    end subroutine vc_initialize


    subroutine vc_update_geometry(self, geometry)
        class(vc_evaluator) :: self
        real(dp), intent(in) :: geometry(:, :)

        if (size(geometry, 1) /= 1) then
            write(stderr, *) 'Error in vc_evaluator_mod, update_geometry subroutine.'
            write(stderr, *) '  Geometry should have dimensions (1, nmode).'
            stop
        end if

        self%q = geometry(1, :)
    end subroutine vc_update_geometry


    subroutine vc_eval(self)
        use linalg_wrapper_mod, only : syev
        class(vc_evaluator) :: self
        type(rmat) :: w_block(3)
        type(rmat) :: w_eigvec(3)
        integer :: i, imode, ns
        real(dp), allocatable :: wrk_e(:)

        ns = self%model%tot_ns

        if (.not. associated(self%model)) then
            write(stderr, *) 'Error in vc_evaluator_mod, eval subroutine.'
            write(stderr, *) '  Model not associated.'
            stop
        end if

        do i = 1, 3
            if (self%model%nstate(i) == 0) cycle
            w_block(i)%c = self%model%w(i)%eval(self%q)
            w_eigvec(i)%c = w_block(i)%c
            if (allocated(wrk_e)) deallocate(wrk_e)
            allocate(wrk_e(size(w_block(i)%c, 1)))
            call syev(w_eigvec(i)%c, wrk_e, jobz='V', uplo='U')
        end do

        if (.not. allocated(self%w_full)) then
             allocate(self%w_full(ns, ns), source=0.0_dp)
        end if
        if (.not. allocated(self%w_eigvec)) then
             allocate(self%w_eigvec(ns, ns), source=0.0_dp)
        end if
        if (.not. allocated(self%w_grad)) then
             allocate(self%w_grad(self%model%nmode, ns, ns), source=0.0_dp)
        end if
        call expand_mult_blocks(w_eigvec, self%w_eigvec)
        call expand_mult_blocks(w_block, self%w_full)
        call change_basis(self%w_eigvec, self%w_full)

        do imode = 1, self%model%nmode
            do i = 1, 3
                if (self%model%nstate(i) == 0) cycle
                w_block(i)%c = self%model%dw(imode, i)%eval(self%q)
            end do
            call expand_mult_blocks(w_block, self%w_grad(imode, :, :))
            call change_basis(self%w_eigvec, self%w_grad(imode, :, :))
        end do

        if (self%model%soc%max_order >= 0) then
            w_block(1)%c = self%model%soc%eval(self%q)
            call change_basis(self%w_eigvec, w_block(1)%c)
            self%w_full = self%w_full + w_block(1)%c
        end if
        if (self%model%dm(1)%max_order >= 0) then
            if (.not. allocated(self%dm)) allocate(self%dm(3, ns, ns))
            self%dm(1, :, :) = self%model%dm(1)%eval(self%q)
            self%dm(2, :, :) = self%model%dm(2)%eval(self%q)
            self%dm(3, :, :) = self%model%dm(3)%eval(self%q)
            call change_basis(self%w_eigvec, self%dm(1, :, :))
            call change_basis(self%w_eigvec, self%dm(2, :, :))
            call change_basis(self%w_eigvec, self%dm(3, :, :))
        end if
    end subroutine vc_eval


    subroutine expand_mult_blocks(blocks, w_full)
        type(rmat), intent(in) :: blocks(:)
        real(dp), intent(out) :: w_full(:, :)
        integer :: i, j, i0, ns

        i0 = 0
        do i = 1, size(blocks)
            if (.not. allocated(blocks(i)%c)) cycle
            ns = size(blocks(i)%c, 1)
            do j = 1, i
                w_full(i0+1:i0+ns, i0+1:i0+ns) = blocks(i)%c
                i0 = i0 + ns
            end do
        end do
    end subroutine expand_mult_blocks

    
    subroutine change_basis(trans, mat)
        real(dp), intent(in) :: trans(:, :)
        real(dp), intent(inout) :: mat(:, :)
        mat = matmul(transpose(trans), mat)
        mat = matmul(mat, trans)
    end subroutine change_basis

    
    function vc_get_nadv(self, istate1, istate2) result(nadv)
        class(vc_evaluator) :: self
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)
        integer :: imode
        real(dp) :: edif
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        allocate(nadv(self%model%nmode), source=0.0_dp)
        edif = self%w_full(istate2, istate2) - self%w_full(istate1, istate1)
        if (abs(edif) < tiny_hf) then
            if (stdp1) then
                write(stderr, *) 'Warning in vc_evaluator_mod, vc_get_nadv function.'
                write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
            end if
            edif = sign(tiny_hf, edif)
        end if
        do imode = 1, self%model%nmode
            nadv(imode) = self%w_grad(imode, istate1, istate2) / edif
        end do
    end function vc_get_nadv

    
    function vc_get_oscill(self, istate) result(oscill)
        class(vc_evaluator) :: self
        integer, intent(in) :: istate
        real(dp), allocatable :: oscill(:)
        integer :: i, jstate
        real(dp) :: edif, mu2

        allocate(oscill(self%model%tot_ns), source=0.0_dp)
        do jstate = 1, self%model%tot_ns
            if (jstate == istate) cycle
            edif = self%w_full(jstate, jstate) - self%w_full(istate, istate)
            mu2 = 0.0_dp
            do i = 1, 3
                if (self%model%dm(i)%max_order >= 0) then
                    mu2 = mu2 + self%dm(i, istate, jstate)**2
                end if
            end do
            oscill(jstate) = 2.0_dp / 3.0_dp * edif * mu2
        end do
    end function vc_get_oscill


end module vc_evaluator_mod