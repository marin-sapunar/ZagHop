module vc_evaluator_mod
    use evaluator_diabatic_mod
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
        procedure :: check_atoms => vc_check_atoms
        procedure :: set_geometry => vc_set_geometry
        procedure :: eval_diab => vc_eval_diab
        procedure :: get_derivative_target_coord => vc_get_derivative_target_coord
        procedure :: get_reference_geometry => vc_get_reference_geometry
        procedure :: default_reference => vc_default_reference
        procedure :: get_oscill => vc_get_oscill
    end type vc_evaluator

contains

    subroutine vc_initialize(self, model)
        class(vc_evaluator) :: self
        type(vc_model), intent(in), target :: model

        self%model => model
        self%classical_coordinates = 'normal_modes'
        self%input_coordinates = self%classical_coordinates
        self%cartesian_frame = 'fixed'
        ! Manifolds of singlets, doublets and triplets. The template lists every component.
        call self%set_spin_manifolds(model%nstate, [1, 2, 3], 'components')
        self%ndof = self%model%nmode
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


    !> Check that the atoms of the system match the reference geometry of the normal modes.
    subroutine vc_check_atoms(self, sym, mass)
        class(vc_evaluator), intent(in) :: self
        character(len=*), intent(in) :: sym(:) !< Atom symbols.
        real(dp), intent(in) :: mass(:) !< Atom masses (amu).

        call self%model%nm%check_atoms(sym, mass)
    end subroutine vc_check_atoms


    !> Normal modes, or a Cartesian geometry that is converted to normal modes. The Cartesian
    !! geometry and the rotation into the Eckart frame are kept for
    !! vc_get_derivative_target_coord.
    subroutine vc_set_geometry(self, geometry, coordinates)
        class(vc_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        if (coordinates == self%classical_coordinates) then
            if (size(geometry, 1) /= 1) then
                write(stderr, *) 'Error in vc_evaluator_mod, set_geometry subroutine.'
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
            write(stderr, *) 'Error in vc_evaluator_mod, set_geometry subroutine.'
            write(stderr, *) '  Unsupported coordinates "'//coordinates//'". Expected "'//      &
            &                self%classical_coordinates//'" or "cartesian".'
            stop
        end if
        self%input_coordinates = coordinates
    end subroutine vc_set_geometry


    !> Convert a derivative with respect to the normal modes to the coordinates of the geometry
    !! passed to set_geometry: unchanged for normal modes, Cartesian for a Cartesian geometry.
    function vc_get_derivative_target_coord(self, d) result(dr)
        class(vc_evaluator), intent(in) :: self
        real(dp), intent(in) :: d(:)
        real(dp), allocatable :: dr(:)

        if (self%input_coordinates /= 'cartesian') then
            dr = d
        else if (self%cartesian_frame == 'eckart') then
            dr = self%model%nm%deriv_nm_to_cart(d, self%cart_geom, self%frame_rot)
        else
            dr = self%model%nm%deriv_nm_to_cart(d)
        end if
    end function vc_get_derivative_target_coord


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


    !> The reference of the template, q = 0.
    subroutine vc_default_reference(self, geom)
        class(vc_evaluator), intent(in) :: self
        real(dp), allocatable, intent(out) :: geom(:, :)

        allocate(geom(1, self%ndof), source=0.0_dp)
    end subroutine vc_default_reference


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

        ! The spin-free blocks are given once for each multiplicity.
        do i = 1, 3
            if (self%model%nstate(i) == 0) cycle
            w_block(i)%c = self%model%w(i)%eval(self%q)
        end do
        self%sf_w = block_diagonal_mat(w_block)

        ! The template only gives constant spin-orbit couplings, so soc_dw stays zero.
        if (self%model%soc%max_order >= 0) then
            self%soc_w = self%model%soc%eval(self%q)
        end if

        do imode = 1, self%model%nmode
            do i = 1, 3
                if (self%model%nstate(i) == 0) cycle
                w_block(i)%c = self%model%dw(imode, i)%eval(self%q)
            end do
            self%sf_dw(imode, :, :) = block_diagonal_mat(w_block)
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
        case default
            write(stderr, *) 'Error in vc_evaluator_mod, vc_get_oscill function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select

        oscill = 2.0_dp / 3.0_dp * edif * sum(dm_bas(:, istate, :)**2, dim=1)
    end function vc_get_oscill


end module vc_evaluator_mod