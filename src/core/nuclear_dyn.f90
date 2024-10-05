!--------------------------------------------------------------------------------------------------
! MODULE: nuclear_dyn_mod
!
! DESCRIPTION: 
!> @brief Subroutines for propagation of nuclear coordinates.
!--------------------------------------------------------------------------------------------------
module nuclear_dyn_mod
    ! Import variables
    use global_defs
    use matrix_mod, only : diagonal_mat
    implicit none

    private
    public :: dyn_updategeom
    public :: dyn_updatevelo
    public :: set_geom_center_of_mass
    public :: project_translation_rotation_from_velocity

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

        masi = diagonal_mat(1.0_dp / mass)
        q = velo - dt * matmul(grad, masi) * 0.5_dp
        if (allocated(ctrl%cns)) call rattle_geom(ctrl%cns, dt, geo1, masi, q)
        geo2 = geo1 + dt * q
        select case(ctrl%orientlvl)
        case(1)
            call set_geom_center_of_mass(mass, geo2)
        case(2)
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

        masi = diagonal_mat(1.0_dp / mass)
        vel2 = vel1 - dt * matmul(grd1 + grd2, masi) * 0.5_dp
        if (ctrl%thermostat == 1) then
            call berendsen_thermostat(ctrl%target_t, ctrl%tau_t, dt, mass, vel2)
        end if
        if (allocated(ctrl%cns)) call rattle_velo(ctrl%cns, geom, masi, vel2)
        select case(ctrl%orientlvl)
        case(2)
            call project_translation_rotation_from_velocity(mass, geom, vel2)
        end select
    end subroutine dyn_updatevelo


    subroutine set_geom_center_of_mass(mass, geom)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: geom(:, :)
        real(dp) :: xyz_cm(size(geom, 1))
        integer :: i

        xyz_cm = matmul(geom, mass) / sum(mass)
        do i = 1, size(mass)
            geom(:, i) = geom(:, i) - xyz_cm
        end do
    end subroutine set_geom_center_of_mass


    subroutine project_translation_rotation_from_velocity(mass, geo, velo)
        !< @todo: Complete cleanup required.
        use matrix_mod, only : mat_inv
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geo(:,:)
        real(dp), intent(inout) :: velo(:,:)
        integer :: i, natom
        real(dp) :: vel_cm(3)
        real(dp) :: ang_mom(3)
        real(dp) :: omega(3)
        real(dp) :: inert_mom(3,3)
        real(dp), allocatable :: velocity(:, :)
        real(dp), allocatable :: geom(:, :)

        velocity = transpose(velo)
        geom = transpose(geo)
        natom = size(mass)
        vel_cm = matmul(mass,velocity) / sum(mass)
        ang_mom = 0.0_dp
        inert_mom = 0.0_dp
        do i = 1, natom
            velocity(i,1:3) = velocity(i,1:3)-vel_cm(1:3)
            ang_mom(1)      = ang_mom(1) + mass(i) * (geom(i,2) * velocity(i,3) - geom(i,3) * velocity(i,2))
            ang_mom(2)      = ang_mom(2) + mass(i) * (geom(i,3) * velocity(i,1) - geom(i,1) * velocity(i,3))
            ang_mom(3)      = ang_mom(3) + mass(i) * (geom(i,1) * velocity(i,2) - geom(i,2) * velocity(i,1))
            inert_mom(1,1)  = inert_mom(1,1) + mass(i) * (geom(i,2)**2 + geom(i,3)**2)
            inert_mom(2,2)  = inert_mom(2,2) + mass(i) * (geom(i,1)**2 + geom(i,3)**2)
            inert_mom(3,3)  = inert_mom(3,3) + mass(i) * (geom(i,1)**2 + geom(i,2)**2)
            inert_mom(1,2)  = inert_mom(1,2) - mass(i) * geom(i,1) * geom(i,2)
            inert_mom(1,3)  = inert_mom(1,3) - mass(i) * geom(i,1) * geom(i,3)
            inert_mom(2,3)  = inert_mom(2,3) - mass(i) * geom(i,2) * geom(i,3)
        end do
        inert_mom(2,1) = inert_mom(1,2)
        inert_mom(3,1) = inert_mom(1,3)
        inert_mom(3,2) = inert_mom(2,3)
        inert_mom = mat_inv(inert_mom)
        omega = matmul(inert_mom, ang_mom)
        do i = 1, natom
            velocity(i,1) = velocity(i,1) - (omega(2) * geom(i,3) - omega(3) * geom(i,2))
            velocity(i,2) = velocity(i,2) - (omega(3) * geom(i,1) - omega(1) * geom(i,3))
            velocity(i,3) = velocity(i,3) - (omega(1) * geom(i,2) - omega(2) * geom(i,1))
        end do
        velo = transpose(velocity)
    end subroutine  project_translation_rotation_from_velocity


    subroutine berendsen_thermostat(target_t, tau_t, dt, mass, vel)
        use system_var, only : ekin
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
