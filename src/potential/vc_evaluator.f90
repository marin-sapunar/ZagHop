module vc_evaluator_mod
    use evaluator_diabatic_mod, diabatic_get_gradient => get_gradient, &
    &                           diabatic_get_nadv => get_nadv
    use vc_model_mod
    use rmat_mod, only : rmat
    implicit none

    type, extends(diabatic_evaluator) :: vc_evaluator
        type(vc_model), pointer :: model => null()
        character(len=:), allocatable :: input_coordinates !< Coordinates of the geometry passed
                                                           !! to update_geometry.
        character(len=:), allocatable :: cartesian_frame !< Frame of Cartesian geometries:
        !! - fixed - Used as given, in the frame of the reference geometry of the normal modes.
        !! - eckart - Rotated into the Eckart frame of the reference geometry.
        real(dp), allocatable :: cart_geom(:, :) !< Cartesian geometry (eckart frame only).
        real(dp) :: frame_rot(3, 3) !< Rotation of cart_geom into the Eckart frame.
        real(dp), allocatable :: q(:)
        
        type(rmat) :: eigvec_block(3)
        type(rmat), allocatable :: grad_block(:, :)
        real(dp), allocatable :: dm(:, :, :)
    contains
        procedure :: initialize => vc_initialize
        procedure :: set_cartesian_frame => vc_set_cartesian_frame
        procedure :: update_geometry => vc_update_geometry
        procedure :: eval_diab => vc_eval_diab
        procedure :: get_gradient => vc_get_gradient
        procedure :: get_nadv => vc_get_nadv
        procedure :: get_reference_geometry => vc_get_reference_geometry
        procedure :: get_oscill => vc_get_oscill
    end type vc_evaluator

contains

    subroutine vc_initialize(self, model)
        class(vc_evaluator) :: self
        type(vc_model), intent(in), target :: model
        integer :: i, i0, i1

        self%model => model
        self%classical_coordinates = 'normal_modes'
        self%input_coordinates = self%classical_coordinates
        self%cartesian_frame = 'fixed'
        ! Blocks of singlets, doublets and triplets.
        call self%set_group_nstate(model%nstate, [1, 2, 3])

        ! The group_adiabatic states are only spin-free if the SOC has no elements inside a group.
        if (model%soc%max_order >= 0) then
            do i = 1, self%n_group
                i0 = self%group_i0(i)
                i1 = i0 + self%group_nstate(i) - 1
                if (any(model%soc%zero(i0:i1, i0:i1) /= 0.0_dp)) self%soc_in_group = .true.
                if (model%soc%max_order >= 1) then
                    if (any(model%soc%linear(:, i0:i1, i0:i1) /= 0.0_dp)) then
                        self%soc_in_group = .true.
                    end if
                end if
                if (model%soc%max_order >= 2) then
                    if (any(model%soc%quad(:, :, i0:i1, i0:i1) /= 0.0_dp)) then
                        self%soc_in_group = .true.
                    end if
                end if
            end do
        end if

        self%ndof = self%model%nmode
        allocate(self%diab_w(self%n_state, self%n_state))
        allocate(self%diab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%adiab_w(self%n_state, self%n_state))
        allocate(self%adiab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%group_adiab_w(self%n_state, self%n_state))
        allocate(self%group_adiab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%group_adiab_trans(self%n_state, self%n_state))
        allocate(self%q(self%ndof))
    end subroutine vc_initialize


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: vc_set_cartesian_frame
    !> @brief Set the frame in which Cartesian geometries are transformed to normal modes.
    !> @details
    !! In the fixed frame, the potential is only invariant to infinitesimal rotations about the
    !! reference geometry, and propagation in Cartesian coordinates is equivalent to propagation
    !! of the normal modes. In the Eckart frame, the potential is invariant to any rotation, so
    !! the angular momentum is conserved and rotation-vibration coupling is included.
    !----------------------------------------------------------------------------------------------
    subroutine vc_set_cartesian_frame(self, frame)
        use rigid_body_mod, only : external_modes
        class(vc_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: frame

        select case(frame)
        case('fixed')
        case('eckart')
            if (allocated(self%model%nm%ref_geom)) then
                if (size(external_modes(self%model%nm%mass, self%model%nm%ref_geom, .true.), 2)  &
                &   /= 6) then
                    write(stderr, *) 'Error in vc_evaluator_mod, set_cartesian_frame subroutine.'
                    write(stderr, *) '  Eckart frame is not implemented for linear molecules.'
                    stop
                end if
            end if
        case default
            write(stderr, *) 'Error in vc_evaluator_mod, set_cartesian_frame subroutine.'
            write(stderr, *) '  Unrecognized frame "'//frame//'". Expected "fixed" or "eckart".'
            stop
        end select
        self%cartesian_frame = frame
    end subroutine vc_set_cartesian_frame


    subroutine vc_update_geometry(self, geometry, coordinates)
        class(vc_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        if (coordinates == self%classical_coordinates) then
            if (size(geometry, 1) /= 1) then
                write(stderr, *) 'Error in vc_evaluator_mod, update_geometry subroutine.'
                write(stderr, *) '  Geometry should have dimensions (1, nmode).'
                stop
            end if
            self%q = geometry(1, :)
        else if (coordinates == 'cartesian') then
            if (self%cartesian_frame == 'eckart') then
                self%cart_geom = geometry
                self%frame_rot = self%model%nm%eckart_rotation(geometry)
                self%q = self%model%nm%cart_to_nm(geometry, self%frame_rot)
            else
                self%q = self%model%nm%cart_to_nm(geometry)
            end if
        else
            write(stderr, *) 'Error in vc_evaluator_mod, update_geometry subroutine.'
            write(stderr, *) '  Unsupported coordinates "'//coordinates//'". Expected "'//      &
            &                self%classical_coordinates//'" or "cartesian".'
            stop
        end if
        self%input_coordinates = coordinates

        if (allocated(self%adiab_trans)) self%ldiab_trans = self%adiab_trans
    end subroutine vc_update_geometry


    !> Gradient in the coordinates of the geometry passed to update_geometry.
    function vc_get_gradient(self, basis, istate) result(grad)
        class(vc_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: grad(:, :)

        grad = diabatic_get_gradient(self, basis, istate)
        if (self%input_coordinates == 'cartesian') then
            grad = reshape(vc_deriv_to_cart(self, grad(1, :)), [3, self%model%nm%natom])
        end if
    end function vc_get_gradient


    !> Nonadiabatic coupling vector in the coordinates of the geometry passed to update_geometry.
    function vc_get_nadv(self, basis, istate1, istate2) result(nadv)
        class(vc_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)

        nadv = diabatic_get_nadv(self, basis, istate1, istate2)
        if (self%input_coordinates == 'cartesian') then
            nadv = vc_deriv_to_cart(self, nadv)
        end if
    end function vc_get_nadv


    !> Convert a derivative with respect to the normal modes to a Cartesian derivative at the
    !! geometry passed to update_geometry.
    function vc_deriv_to_cart(self, dq) result(dr)
        class(vc_evaluator), intent(in) :: self
        real(dp), intent(in) :: dq(:)
        real(dp), allocatable :: dr(:)

        if (self%cartesian_frame == 'eckart') then
            dr = self%model%nm%deriv_nm_to_cart(dq, self%cart_geom, self%frame_rot)
        else
            dr = self%model%nm%deriv_nm_to_cart(dq)
        end if
    end function vc_deriv_to_cart


    !> In Cartesian coordinates in the fixed frame, the normal modes are defined relative to the
    !! reference geometry from the V0 file. In the Eckart frame, the potential is invariant to
    !! rotation.
    subroutine vc_get_reference_geometry(self, coordinates, geom)
        class(vc_evaluator), intent(in) :: self
        character(len=*), intent(in) :: coordinates
        real(dp), allocatable, intent(out) :: geom(:, :)

        if (coordinates /= 'cartesian' .or. self%cartesian_frame /= 'fixed') return
        if (allocated(self%model%nm%ref_geom)) geom = self%model%nm%ref_geom
    end subroutine vc_get_reference_geometry


    subroutine vc_eval_diab(self)
        use matrix_mod, only : block_diagonal_mat
        class(vc_evaluator), intent(inout) :: self
        type(rmat) :: w_block(3)
        integer :: i, imode

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

        do imode = 1, self%model%nmode
            do i = 1, 3
                if (self%model%nstate(i) == 0) cycle
                w_block(i)%c = self%model%dw(imode, i)%eval(self%q)
            end do
            self%diab_dw(imode, :, :) = block_diagonal_mat(w_block, [1, 2, 3])
        end do

        if (self%model%dm(1)%max_order >= 0) then
            if (.not. allocated(self%dm)) allocate(self%dm(3, self%n_state, self%n_state))
            self%dm(1, :, :) = self%model%dm(1)%eval(self%q)
            self%dm(2, :, :) = self%model%dm(2)%eval(self%q)
            self%dm(3, :, :) = self%model%dm(3)%eval(self%q)
        end if
    end subroutine vc_eval_diab

    
    function vc_get_oscill(self, basis, istate) result(oscill)
        class(vc_evaluator) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: oscill(:)
        integer :: jstate
        real(dp), allocatable :: edif(:)
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

        oscill = 2.0_dp / 3.0_dp * edif * sum(dm_bas(:, istate, :)**2, dim=1)
    end function vc_get_oscill


end module vc_evaluator_mod