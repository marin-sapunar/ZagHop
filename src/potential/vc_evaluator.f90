module vc_evaluator_mod
    use evaluator_diabatic_mod
    use vc_model_mod
    use rmat_mod, only : rmat
    implicit none

    type, extends(diabatic_evaluator) :: vc_evaluator
        type(vc_model), pointer :: model => null()
        real(dp), allocatable :: q(:)
        
        type(rmat) :: eigvec_block(3)
        type(rmat), allocatable :: grad_block(:, :)
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
        integer :: i, cgroup, ci, j

        self%model => model
        self%n_group = sum([1, 2, 3], mask=(model%nstate > 0))
        allocate(self%group_nstate(self%n_group))
        allocate(self%group_i0(self%n_group))
        cgroup = 1
        ci = 1
        do i = 1, 3
            if (model%nstate(i) == 0) cycle
            do j = 1, i
                self%group_nstate(cgroup) = model%nstate(i)
                self%group_i0(cgroup) = ci
                cgroup = cgroup + 1
                ci = ci + model%nstate(i)
            end do
        end do
        self%n_state = sum(self%group_nstate)
        allocate(self%diab_w(self%n_state, self%n_state))
        allocate(self%diab_dw(self%model%nmode, self%n_state, self%n_state))
        allocate(self%adiab_w(self%n_state, self%n_state))
        allocate(self%adiab_dw(self%model%nmode, self%n_state, self%n_state))
        allocate(self%adiab_trans(self%n_state, self%n_state))
        allocate(self%group_adiab_w(self%n_state, self%n_state))
        allocate(self%group_adiab_dw(self%model%nmode, self%n_state, self%n_state))
        allocate(self%group_adiab_trans(self%n_state, self%n_state))
        allocate(self%q(self%model%nmode))
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


    subroutine vc_eval(self, previous)
        use linalg_wrapper_mod, only : syev
        use matrix_mod, only : block_diagonal_mat
        class(vc_evaluator) :: self
        real(dp), intent(in), optional :: previous(:, :)
        type(rmat) :: w_block(3)
        integer :: i, imode
        real(dp), allocatable :: wrk_mat(:, :)
        integer, allocatable :: phase(:)
        real(dp), allocatable :: eigvec_olap(:, :)

        if (.not. associated(self%model)) then
            write(stderr, *) 'Error in vc_evaluator_mod, eval subroutine.'
            write(stderr, *) '  Model not associated.'
            stop
        end if

        do i = 1, 3
            if (self%model%nstate(i) == 0) cycle
            w_block(i)%c = self%model%w(i)%eval(self%q)
        end do
        self%diab_w = block_diagonal_mat(w_block, [1, 2, 3])

        if (self%model%soc%max_order >= 0) then
            self%diab_w = self%diab_w + self%model%soc%eval(self%q)
        end if

        call self%eval_eigvec('adiabatic')
        call self%eval_eigvec('group_adiabatic')

        do imode = 1, self%model%nmode
            do i = 1, 3
                if (self%model%nstate(i) == 0) cycle
                w_block(i)%c = self%model%dw(imode, i)%eval(self%q)
            end do
            self%diab_dw(imode, :, :) = block_diagonal_mat(w_block, [1, 2, 3])
        end do
        
        if (present(previous)) then
            call match_phase(previous, self%adiab_trans)
          !  call match_phase(previous%group_adiab_trans, self%group_adiab_trans)
        end if

        self%adiab_w = self%diab_w
        self%adiab_dw = self%diab_dw
        call self%change_basis('adiabatic', self%adiab_w)
        do imode = 1, self%model%nmode
            call self%change_basis('adiabatic', self%adiab_dw(imode, :, :))
        end do

        self%group_adiab_w = self%diab_w
        self%group_adiab_dw = self%diab_dw
        call self%change_basis('group_adiabatic', self%group_adiab_w)
        do imode = 1, self%model%nmode
            call self%change_basis('group_adiabatic', self%group_adiab_dw(imode, :, :))
        end do

        if (self%model%dm(1)%max_order >= 0) then
            if (.not. allocated(self%dm)) allocate(self%dm(3, self%n_state, self%n_state))
            self%dm(1, :, :) = self%model%dm(1)%eval(self%q)
            self%dm(2, :, :) = self%model%dm(2)%eval(self%q)
            self%dm(3, :, :) = self%model%dm(3)%eval(self%q)
        end if
    end subroutine vc_eval

    
    function vc_get_nadv(self, basis, istate1, istate2) result(nadv)
        class(vc_evaluator) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)
        integer :: imode
        real(dp) :: edif
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        allocate(nadv(self%model%nmode), source=0.0_dp)
        select case(basis)
        case('diabatic')
            continue
        case('adiabatic')
            edif = self%adiab_w(istate2, istate2) - self%adiab_w(istate1, istate1)
            if (abs(edif) < tiny_hf) then
                if (stdp1) then
                    write(stderr, *) 'Warning in vc_evaluator_mod, vc_get_nadv function.'
                    write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                    write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
                end if
                edif = sign(tiny_hf, edif)
            end if
            do imode = 1, self%model%nmode
                nadv(imode) = self%adiab_dw(imode, istate1, istate2) / edif
            end do
        case('group_adiabatic')
             edif = self%group_adiab_w(istate2, istate2) - self%group_adiab_w(istate1, istate1)
             if (abs(edif) < tiny_hf) then
                if (stdp1) then
                    write(stderr, *) 'Warning in vc_evaluator_mod, vc_get_nadv function.'
                    write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                    write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
                end if
                edif = sign(tiny_hf, edif)
            end if
            do imode = 1, self%model%nmode
                nadv(imode) = self%group_adiab_dw(imode, istate1, istate2) / edif
            end do
        case default
            write(stderr, *) 'Error in vc_evaluator_mod, vc_get_nadv function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function vc_get_nadv

    
    function vc_get_oscill(self, basis, istate) result(oscill)
        class(vc_evaluator) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: oscill(:)
        integer :: i, jstate
        real(dp), allocatable :: edif(:)
        real(dp) :: mu2
        real(dp), allocatable :: dm_bas(:, :, :)

        allocate(oscill(self%model%tot_ns), source=0.0_dp)
        allocate(edif(self%model%tot_ns), source=0.0_dp)
        dm_bas = self%dm
        call self%change_basis(basis, dm_bas(1, :, :))
        call self%change_basis(basis, dm_bas(2, :, :))
        call self%change_basis(basis, dm_bas(3, :, :))
        select case(basis)
        case('diabatic')
            do jstate = 1, self%model%tot_ns
                edif(jstate) = self%diab_w(jstate, jstate) - self%diab_w(istate, istate)
            end do
        case('adiabatic')
            do jstate = 1, self%model%tot_ns
                edif(jstate) = self%adiab_w(jstate, jstate) - self%adiab_w(istate, istate)
            end do
        case('group_adiabatic')
            do jstate = 1, self%model%tot_ns
                edif(jstate) = self%group_adiab_w(jstate, jstate) - self%group_adiab_w(istate, istate)
            end do
        case default
            write(stderr, *) 'Error in vc_evaluator_mod, vc_get_oscill function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select

        oscill = 2.0_dp / 3.0_dp * edif * sum(dm_bas(:, istate, :)**2, dim=2)
    end function vc_get_oscill


end module vc_evaluator_mod