!--------------------------------------------------------------------------------------------------
! MODULE: nuclear_dyn_mod
!
! DESCRIPTION: 
!> @brief Subroutines for propagation of nuclear coordinates.
!--------------------------------------------------------------------------------------------------
module nuclear_dyn_mod
    ! Import variables
    use global_defs
    use matrix_mod, only : diag
    implicit none

    private
    public :: dyn_updategeom
    public :: dyn_updatevelo
    public :: dyn_initial_orientation
    public :: set_geom_center_of_mass
    public :: project_external_velocity

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: dyn_updategeom
    !
    ! DESCRIPTION:
    !> @brief Update geometry using the velocity Verlet algorithm.
    !> @details
    !! When updating the geometry, the rattle algorithm is called to ensure any constraints are
    !! satisfied.
    !----------------------------------------------------------------------------------------------
    subroutine dyn_updategeom(dt, mass, geo1, grad, velo, geo2)
        use control_var
        use rattle_mod
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geo1(:, :)
        real(dp), intent(in) :: grad(:, :)
        real(dp), intent(in) :: velo(:, :)
        real(dp), intent(out) :: geo2(:, :)
        real(dp) :: masi(size(mass), size(mass))
        real(dp) :: q(size(geo2, 1), size(geo1, 2))

        if (stdp3) write(stdout, *) '  Updating geometry using velocity Verlet algorithm.'

        masi = diag(1.0_dp / mass)
        q = velo - dt * matmul(grad, masi) * 0.5_dp
        if (allocated(ctrl%cns)) call rattle_geom(ctrl%cns, dt, geo1, masi, q)
        geo2 = geo1 + dt * q
        select case(ctrl%orientlvl)
        case(1, 2)
            call set_geom_center_of_mass(mass, geo2)
        end select
    end subroutine dyn_updategeom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: dyn_updatevelo
    !
    ! DESCRIPTION:
    !> @brief Update velocity using the velocity Verlet algorithm.
    !> @details
    !! After updating the velocity, the rattle algorithm is called to ensure any constraints are
    !! satisfied.
    !----------------------------------------------------------------------------------------------
    subroutine dyn_updatevelo(dt, mass, geom, grd1, grd2, vel1, vel2)
        use control_var
        use rattle_mod
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :)
        real(dp), intent(in) :: grd1(:, :)
        real(dp), intent(in) :: grd2(:, :)
        real(dp), intent(in) :: vel1(:, :)
        real(dp), intent(out) :: vel2(:, :)
        real(dp) :: masi(size(mass), size(mass))

        if (stdp3) write(stdout, *) '  Updating velocity using velocity Verlet algorithm.'

        masi = diag(1.0_dp / mass)
        vel2 = vel1 - dt * matmul(grd1 + grd2, masi) * 0.5_dp
        if (ctrl%thermostat == 1) then
            call berendsen_thermostat(ctrl%target_t, ctrl%tau_t, dt, mass, vel2)
        end if
        if (allocated(ctrl%cns)) call rattle_velo(ctrl%cns, geom, masi, vel2)
        call remove_external_motion(mass, geom, vel2)
    end subroutine dyn_updatevelo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: dyn_initial_orientation
    !
    ! DESCRIPTION:
    !> @brief Remove overall translation and rotation from the initial conditions.
    !> @details
    !! The geometry is moved to the centre of mass and translation (orientlvl 1) or translation
    !! and rotation (orientlvl 2) are projected from the velocity. For orientlvl 2, if the
    !! potential is defined relative to a reference geometry (ctrl%orient_ref), the geometry and
    !! velocity are first rotated into the frame of the reference geometry (Eckart frame).
    !----------------------------------------------------------------------------------------------
    subroutine dyn_initial_orientation(mass, geom, velo)
        use control_var
        use rigid_body_mod, only : kabsch_rotation
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: geom(:, :)
        real(dp), intent(inout) :: velo(:, :)
        real(dp) :: u(3, 3)
        real(dp) :: angle

        if (ctrl%orientlvl == 0) return
        call set_geom_center_of_mass(mass, geom)
        if (ctrl%orientlvl == 2 .and. allocated(ctrl%orient_ref)) then
            u = kabsch_rotation(mass, geom, ctrl%orient_ref)
            geom = matmul(u, geom)
            velo = matmul(u, velo)
            if (stdp2) then
                angle = acos(max(-1.0_dp, min(1.0_dp, 0.5_dp * (u(1,1) + u(2,2) + u(3,3) - 1))))
                write(stdout, '(1x,a,f0.4,a)') 'Initial geometry rotated by ', angle * 180 / pi,  &
                &                              ' degrees into the frame of the reference geometry.'
            end if
        end if
        call remove_external_motion(mass, geom, velo)
    end subroutine dyn_initial_orientation


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: remove_external_motion
    !
    ! DESCRIPTION:
    !> @brief Project translation (orientlvl 1) or translation and rotation (orientlvl 2) from
    !! the velocity.
    !> @details
    !! If the potential is defined relative to a reference geometry (ctrl%orient_ref), rotations
    !! are taken about the reference geometry. Such a potential is only invariant to infinitesimal
    !! rotations about the reference geometry, so these (and not rotations about the current
    !! geometry) are the motions which do not couple to the rest of the dynamics.
    !----------------------------------------------------------------------------------------------
    subroutine remove_external_motion(mass, geom, velo)
        use control_var
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :)
        real(dp), intent(inout) :: velo(:, :)

        select case(ctrl%orientlvl)
        case(1)
            call project_external_velocity(mass, geom, .false., velo)
        case(2)
            if (allocated(ctrl%orient_ref)) then
                call project_external_velocity(mass, ctrl%orient_ref, .true., velo)
            else
                call project_external_velocity(mass, geom, .true., velo)
            end if
        end select
    end subroutine remove_external_motion


    subroutine set_geom_center_of_mass(mass, geom)
        use rigid_body_mod, only : center_of_mass
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: geom(:, :)

        geom = geom - spread(center_of_mass(mass, geom), 2, size(geom, 2))
    end subroutine set_geom_center_of_mass


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: project_external_velocity
    !
    ! DESCRIPTION:
    !> @brief Remove overall translation and, optionally, rotation from the velocity.
    !> @details
    !! The mass-weighted velocity is projected onto the orthogonal complement of the translation
    !! and rotation vectors, with the rotations taken about the centre of mass of geom. This
    !! removes the linear and angular momentum with respect to geom while keeping the rest of the
    !! velocity unchanged.
    !----------------------------------------------------------------------------------------------
    subroutine project_external_velocity(mass, geom, rotation, velo)
        use rigid_body_mod, only : external_modes
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :) !< Geometry defining the rotation vectors.
        logical, intent(in) :: rotation !< Also project rotation.
        real(dp), intent(inout) :: velo(:, :)
        real(dp), allocatable :: ext(:, :)
        real(dp) :: sqrt_mass(size(velo))
        real(dp) :: x(size(velo))

        ext = external_modes(mass, geom, rotation)
        sqrt_mass = reshape(spread(sqrt(mass), 1, size(velo, 1)), [size(velo)])
        x = reshape(velo, [size(velo)]) * sqrt_mass
        x = x - matmul(ext, matmul(x, ext))
        velo = reshape(x / sqrt_mass, shape(velo))
    end subroutine project_external_velocity


    subroutine berendsen_thermostat(target_t, tau_t, dt, mass, vel)
        use system_type_mod, only : ekin
        real(dp), intent(in) :: tau_t
        real(dp), intent(in) :: target_t
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: vel(:, :)
        real(dp):: bt_lambda
        real(dp):: current_t

        current_t = 315773.0_dp * f2o3 * ekin(mass, vel) / size(mass)
        bt_lambda = sqrt(num1 + (dt / tau_t) * ((target_t / current_t) - num1))
        vel = vel * bt_lambda
    end subroutine berendsen_thermostat


end module nuclear_dyn_mod
