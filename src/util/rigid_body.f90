!--------------------------------------------------------------------------------------------------
! MODULE: rigid_body_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2026
!
! DESCRIPTION:
!> @brief Overall translation and rotation of a molecule.
!> @details
!! Geometries have shape (3, natom). The masses only enter as relative weights, so any mass unit
!! can be used.
!--------------------------------------------------------------------------------------------------
module rigid_body_mod
    use global_defs
    implicit none

    private
    public :: center_of_mass
    public :: external_modes
    public :: kabsch_rotation


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: center_of_mass
    !----------------------------------------------------------------------------------------------
    pure function center_of_mass(mass, geom) result(com)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :)
        real(dp) :: com(size(geom, 1))

        com = matmul(geom, mass) / sum(mass)
    end function center_of_mass


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: external_modes
    !> @brief Orthonormal translation and rotation vectors in mass-weighted coordinates.
    !> @details
    !! Returns the vectors spanning overall translation and, optionally, infinitesimal rotation
    !! about the centre of mass of geom, in mass-weighted coordinates x = M^(1/2) R. The vectors
    !! are orthonormalized with the modified Gram-Schmidt procedure. Rotation vectors which are
    !! (numerically) linearly dependent on the previous vectors are dropped, so only 5 vectors
    !! are returned for linear molecules and 3 for a single atom.
    !----------------------------------------------------------------------------------------------
    function external_modes(mass, geom, rotation) result(ext)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :) !< Geometry, shape (3, natom).
        logical, intent(in) :: rotation !< Include rotation vectors.
        real(dp), allocatable :: ext(:, :) !< Orthonormal vectors, shape (3*natom, next).
        real(dp), parameter :: tol = 1.0e-6_dp !< Relative norm below which a vector is dropped.
        real(dp) :: sqrt_mass(size(mass))
        real(dp) :: r(3, size(mass))
        real(dp) :: wrk(3, size(mass))
        real(dp) :: vec(3 * size(mass), 6)
        real(dp) :: ref_norm
        integer :: natom, nvec, next, i, j, k, pass

        natom = size(mass)
        sqrt_mass = sqrt(mass)
        r = geom - spread(center_of_mass(mass, geom), 2, natom)
        nvec = 3
        if (rotation) nvec = 6

        ! Translation along axis k.
        do k = 1, 3
            wrk = 0.0_dp
            wrk(k, :) = sqrt_mass
            vec(:, k) = reshape(wrk, [3 * natom])
        end do
        ! Rotation about axis k: e_k x r.
        if (rotation) then
            do k = 1, 3
                i = mod(k, 3) + 1
                j = mod(k + 1, 3) + 1
                wrk = 0.0_dp
                wrk(j, :) = r(i, :) * sqrt_mass
                wrk(i, :) = -r(j, :) * sqrt_mass
                vec(:, 3 + k) = reshape(wrk, [3 * natom])
            end do
        end if

        ref_norm = maxval(norm2(vec(:, 1:nvec), dim=1))
        next = 0
        do k = 1, nvec
            ! Two passes of modified Gram-Schmidt for numerical stability.
            do pass = 1, 2
                do j = 1, next
                    vec(:, k) = vec(:, k) - dot_product(vec(:, j), vec(:, k)) * vec(:, j)
                end do
            end do
            if (norm2(vec(:, k)) < tol * ref_norm) cycle
            next = next + 1
            vec(:, next) = vec(:, k) / norm2(vec(:, k))
        end do
        ext = vec(:, 1:next)
    end function external_modes


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: kabsch_rotation
    !> @brief Rotation which best aligns a geometry with a reference geometry.
    !> @details
    !! Returns the proper rotation matrix U which minimizes sum_i m_i |U r_i - r0_i|^2, where r_i
    !! and r0_i are the positions relative to the centres of mass of geom and ref. At the minimum,
    !! sum_i m_i (U r_i) x r0_i = 0, which are the Eckart rotational conditions with ref as the
    !! reference geometry.
    !----------------------------------------------------------------------------------------------
    function kabsch_rotation(mass, geom, ref) result(u)
        use linalg_wrapper_mod, only : gesvd
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: geom(:, :) !< Geometry to rotate, shape (3, natom).
        real(dp), intent(in) :: ref(:, :) !< Reference geometry, shape (3, natom).
        real(dp) :: u(3, 3)
        real(dp) :: r(3, size(mass))
        real(dp) :: r0(3, size(mass))
        real(dp) :: h(3, 3), w(3, 3), vt(3, 3), s(3)
        real(dp) :: d(3, 3)
        integer :: natom

        natom = size(mass)
        r = geom - spread(center_of_mass(mass, geom), 2, natom)
        r0 = ref - spread(center_of_mass(mass, ref), 2, natom)

        ! Mass-weighted covariance matrix h = sum_i m_i r_i r0_i^T = W S V^T.
        h = matmul(r * spread(mass, 1, 3), transpose(r0))
        call gesvd(h, s, u=w, vt=vt)

        ! U = V D W^T, where D = diag(1, 1, det(V W^T)) excludes reflections.
        d = 0.0_dp
        d(1, 1) = 1.0_dp
        d(2, 2) = 1.0_dp
        d(3, 3) = sign(1.0_dp, det3(vt) * det3(w))
        u = matmul(transpose(vt), matmul(d, transpose(w)))
    end function kabsch_rotation


    pure function det3(a) result(det)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: det

        det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))                                   &
        &   - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))                                   &
        &   + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
    end function det3


end module rigid_body_mod
