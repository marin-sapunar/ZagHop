!--------------------------------------------------------------------------------------------------
! MODULE: rattle_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
! DESCRIPTION: 
!> @brief RATTLE algorithm for constraints in dynamics calculations.
!> @details
!! Subroutines for updating coordinates and velocities for dynamics with constrained bond lengths!! using the RATTLE algorithm.
!! References:
!! - H. C. Andersen: J. Comput. Phys. 52 (1983) 24-34
!! - - DOI: 10.1016/0021-9991(83)90014-1
!--------------------------------------------------------------------------------------------------

module rattle_mod
    use global_defs
    implicit none


    type constraint
        integer :: i(4) !< Indexes of constrained atoms.
        real(dp) :: cval !< Constrained value.
        real(dp) :: tol !< Tolerance.
    end type constraint
    

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Rattle_Geom
    !
    ! DESCRIPTION:
    !> @brief Update displacements for new geometry to match constraints.
    !----------------------------------------------------------------------------------------------
    subroutine rattle_geom(r, dt, geom, im, q)
        type(constraint), intent(in) :: r(:) !< Constraints.
        real(dp), intent(in) :: dt !< Time step.
        real(dp), intent(in) :: geom(:, :) !< Old geometry.
        real(dp), intent(in) :: im(:, :) !< Inverted mass matrix.
        real(dp), intent(inout) :: q(:, :) !< Current displacements.
        real(dp) :: rij(size(q, 1)) !< Old bond vector.
        real(dp) :: s(size(q, 1)) !< Bond vector with displacements.
        real(dp) :: cons !< Constraint test.
        real(dp) :: g !< Change in displacement.
        integer :: c !< Current constraint.
        integer :: i !< Index of bond atom 1.
        integer :: j !< Index of bond atom 2.
        integer :: iter !< Number of RATTLE iterations.
        integer, parameter :: maxiter = 100 !< Fail after # iterations.

        iter = 0
 100    c = 0 
 101    c = c + 1
        i = r(c)%i(1)
        j = r(c)%i(2)
        rij = geom(:, i) - geom(:, j)
        s = rij + dt * (q(:, i) - q(:, j))
        cons = dot_product(s,s) - r(c)%cval**2
        if (abs(cons) < r(c)%tol) then
            if (c == size(r)) return
            goto 101
        end if
        g = cons / (im(i, i) + im(j, j)) / 2 / dt / dot_product(s, rij)
        q(:, i) = q(:, i) - g * rij * im(i, i)
        q(:, j) = q(:, j) + g * rij * im(j, j)

        iter = iter + 1
        if (iter > maxiter) then
            write(stderr, *) 'Error in RattleMod. Maximum number of iterations exceeded.'
            stop
        end if
        goto 100
    end subroutine rattle_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Rattle_Velo
    !
    ! DESCRIPTION:
    !> @brief Update velocities to match constraints.
    !----------------------------------------------------------------------------------------------
    subroutine rattle_velo(r, geom, im, velo)
        type(constraint), intent(in) :: r(:) !< Constraints.
        real(dp), intent(in) :: geom(:, :) !< Current geometry.
        real(dp), intent(in) :: im(:, :) !< Inverted mass matrix.
        real(dp), intent(inout) :: velo(:, :) !< Current velocity.
        real(dp) :: rij(size(velo, 1)) !< Current bond vector.
        real(dp) :: vij(size(velo, 1)) !< Velocity along bond vector
        real(dp) :: cons !< Constraint test.
        real(dp) :: k !< Change in velocities.
        integer :: c !< Current constraint.
        integer :: i !< Index of bond atom 1.
        integer :: j !< Index of bond atom 2.
        integer :: iter !< Number of RATTLE iterations.
        integer, parameter :: maxiter = 100 !< Fail after # iterations.

        iter = 0
 100    c = 0 
 101    c = c + 1
        i = r(c)%i(1)
        j = r(c)%i(2)
        rij = geom(:, i) - geom(:, j)
        vij = velo(:, i) - velo(:, j)
        cons = dot_product(rij, vij)
        if (abs(cons) < r(c)%tol) then
            if (c == size(r)) return
            goto 101
        end if
        k = cons / r(c)%cval**2 / (im(i, i) + im(j, j))
        velo(:, i) = velo(:, i) - k * rij * im(i, i)
        velo(:, j) = velo(:, j) + k * rij * im(j, j)

        iter = iter + 1
        if (iter > maxiter) then
            write(stderr, *) 'Error in RattleMod. Maximum number of velocity iterations exceeded.'
            stop
        end if
        goto 100
    end subroutine rattle_velo


end module rattle_mod
