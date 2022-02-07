!--------------------------------------------------------------------------------------------------
! MODULE: orthog_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
!> @brief Simple orthogonalization subroutines.
!--------------------------------------------------------------------------------------------------
module orthog_mod
    use global_defs
    implicit none


    private
    public :: orthog_gs
    public :: orthog_lowdin


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: orthog_gs
    !
    !> @brief Modified Gram-Schmidt orthogonalization of a set of column vectors.
    !> @details
    !! Input is a real matrix with dimensions n x m matrix containing m vectors of dimension n.
    !! The vectors are orthogonalized (in place) using the Gram-Schmidt algorithm.
    !----------------------------------------------------------------------------------------------
    subroutine orthog_gs(u)
        real(dp), intent(inout) :: u(:, :) !< Column vectors.
        integer :: j
        integer :: i

        u(:, 1) = u(:, 1) / norm(u(:, 1))
        do i = 2, size(u, 2)
            do j = 1, i - 1
                u(:, i) = u(:, i) - proj(u(:, i), u(:, j)) * u(:, j)
            end do
            u(:, i) = u(:, i) / norm(u(:, i))
        end do
    end subroutine orthog_gs

    pure function proj(v1, v2) result(p)
        real(dp), intent(in) :: v1(:)
        real(dp), intent(in) :: v2(:)
        real(dp) :: p
        p = dot_product(v1, v2) / dot_product(v2, v2)
    end function proj

    pure function norm(v) result(n)
        real(dp), intent(in) :: v(:)
        real(dp) :: n
        n = sqrt(dot_product(v, v))
    end function norm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: orthog_lowdin
    !
    !> @brief Lowdin (symmetric) orthogonalization of a set of column vectors.
    !> @details
    !! Input is a real matrix with dimensions n x m matrix containing m vectors of dimension n.
    !! The vectors are orthogonalized (in place) using singular value decomposition.
    !----------------------------------------------------------------------------------------------
    subroutine orthog_lowdin(mat)
        use linalg_wrapper_mod, only : gesvd, gemm
        real(dp), intent(inout) :: mat(:, :) !< Column vectors.
        real(dp), allocatable :: s(:) !< Singular values (not used.)
        real(dp), allocatable :: vt(:, :) !< Transposed V matrix.
        real(dp), allocatable :: u(:, :) !< U matrix.
        integer :: m, n !< Array size.
        integer :: info !< LAPACK return code.

        m = size(mat, 1)
        n = size(mat, 2)
        allocate(s(min(m, n)))
        allocate(u(m, min(m, n)))
        allocate(vt(n, min(m, n)))
        
        ! Compute SVD
        call gesvd(mat, s, u=u, vt=vt, info=info)
        if (info /= 0) call errstop('ORTH_LOWDIN', 'GESVD call failed.', info)

        ! Overwrite input matrix with orthogonized matrix U.Vt.
        call gemm(u, vt, mat)
    end subroutine orthog_lowdin


end module orthog_mod
