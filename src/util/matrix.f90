!--------------------------------------------------------------------------------------------------
! MODULE: matrix_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief LAPACK interfaces and common matrix operations.
!--------------------------------------------------------------------------------------------------
module matrix_mod
    use global_defs, only : dp, errstop
    use linalg_wrapper_mod
    implicit none


    private
    public :: diagonal_mat
    public :: unit_mat
    public :: vec_outer
    public :: mat_norm
    public :: mat_inv
    public :: mat_sy_ev
    public :: mat_sy_exp
    public :: mat_ge_det
    public :: mat_ge_mmm


    interface vec_outer
        module procedure vec_outer_int
        module procedure vec_outer_dp
    end interface vec_outer


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: diagonal_mat
    !> @brief Given a vector return diagonal matrix with values from vector on the diagonal.
    !----------------------------------------------------------------------------------------------
    function diagonal_mat(vec) result(mat)
        real(dp), intent(in) :: vec(:)
        real(dp) :: mat(size(vec), size(vec))
        integer :: i

        mat = num0
        do i = 1, size(vec)
            mat(i, i) = vec(i)
        end do
    end function diagonal_mat


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: unit_mat
    !> @brief Return a unit matrix of dimension n.
    !----------------------------------------------------------------------------------------------
    function unit_mat(n) result(mat)
        integer, intent(in) :: n
        real(dp) :: mat(n, n)
        integer :: i

        mat = num0
        do i = 1, n
            mat(i, i) = num1
        end do
    end function unit_mat


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: vec_outer_int
    !> @brief Give outer product of two int vectors.
    !----------------------------------------------------------------------------------------------
    function vec_outer_int(vec1, vec2) result(mat)
        integer :: vec1(:)
        integer :: vec2(:)
        integer :: mat(size(vec1), size(vec2))
        integer :: i

        do i = 1, size(vec2)
            mat(:, i) = vec2(i) * vec1
        end do
    end function vec_outer_int


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: vec_outer_dp
    !> @brief Give outer product of two dp vectors.
    !----------------------------------------------------------------------------------------------
    function vec_outer_dp(vec1, vec2) result(mat)
        real(dp) :: vec1(:)
        real(dp) :: vec2(:)
        real(dp) :: mat(size(vec1), size(vec2))

        mat = 0.0_dp
        call ger(mat, vec1, vec2)
    end function vec_outer_dp


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_norm
    !> @brief Returns norm of matrix (Frobenius or spectral).
    !----------------------------------------------------------------------------------------------
    function mat_norm(mat, norm_type) result(norm)
        real(dp), intent(in) :: mat(:, :) !< Matrix.
        character(len=*), intent(in), optional :: norm_type !< Type of norm requested.
        real(dp) :: norm !< Norm of matrix.
        real(dp), allocatable :: s(:) !< Singular values.
        real(dp), allocatable :: a(:, :) !< Temporary matrix.
        character(len=:), allocatable :: opt

        opt = 'frobenius'
        if (present(norm_type)) opt = norm_type

        select case(opt)
        case('frobenius')
            norm = sqrt(sum(mat**2))
        case('spectral')
            allocate(s(min(size(mat, 1), size(mat, 2))))
            allocate(a, source=mat)
            call gesvd(a, s)
            norm = sqrt(s(1))
        end select
    end function mat_norm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_inv
    !> @brief Invert of a general matrix
    !----------------------------------------------------------------------------------------------
    function mat_inv(a) result(a_inv)
        real(dp), intent(in) :: a(:, :) !< Input matrix.
        real(dp), allocatable :: a_inv(:, :)
        integer :: ipiv(size(a, 2))
        integer :: n
        integer :: info

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        a_inv = a

        call getrf(a_inv, ipiv=ipiv, info=info)
        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        call getri(a_inv, ipiv, info=info)
        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
    end function mat_inv


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_sy_ev
    !> @brief Eigensystem of a symmetric real matrix.
    !----------------------------------------------------------------------------------------------
    subroutine mat_sy_ev(a, evec, eval)
        real(dp), intent(in) :: a(:, :) !< Input matrix.
        real(dp), allocatable, intent(out) :: evec(:, :) !< Eigenvectors of the symmetric matrix.
        real(dp), allocatable, intent(out) :: eval(:) !< Eigenvalues of the symmetric matrix.
        integer :: n
        integer :: info

        n = size(a, 1)
        if (allocated(evec)) deallocate(evec)
        allocate(evec, source=a)
        if (.not. allocated(eval)) allocate(eval(n))
 
        ! Calculate eigensystem.
        call syev(a=evec, w=eval, jobz='V', info=info)
        if (info /= 0) call errstop('mat_sy_ev', 'syev call failed.', info)
    end subroutine mat_sy_ev


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: mat_sy_exp
    !
    !> @brief Matrix exponential of a symmetric real matrix multiplied by a complex factor.
    !> @details
    !! The subroutine first calculates the eigenvectors/eigenvalues of the matrix. Then the matrix
    !! exponantial is calculated in the eigenvector basis and the matrix is converted back to the
    !! original basis.
    !----------------------------------------------------------------------------------------------
    function mat_sy_exp(a, mult) result (res)
        real(dp), intent(in) :: a(:, :) !< Input matrix.
        complex(dp), intent(in) :: mult !< Multiplier in exponent.
        complex(dp), allocatable :: res(:, :) !< Result matrix.
        real(dp), allocatable :: evec(:, :) !< Eigenvectors of the symmetric matrix.
        real(dp), allocatable :: eval(:) !< Eigenvalues of the symmetric matrix.
        complex(dp), dimension(:, :), allocatable :: tmp1, tmp2 !< Work matrices
        integer :: i, n
        integer :: info

        n = size(a, 1)
        allocate(evec, source=a)
        allocate(eval(n))
        allocate(tmp1(n, n))
        allocate(tmp2(n, n))
        if (.not. allocated(res)) allocate(res(n, n))
 
        ! Calculate eigensystem.
        call syev(a=evec, w=eval, jobz='V', info=info)
        if (info /= 0) call errstop('mat_sy_exp', 'syev call failed.', info)

        ! Matrix exponential in eigenvector basis.
        res = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            res(i, i) = exp(eval(i) * mult)
        end do

        ! Conversion to original basis.
        tmp1 = cmplx(evec, 0.0_dp, dp)
        call gemm(tmp1, res, tmp2)
        call gemm(tmp2, tmp1, res, transb='C')
    end function mat_sy_exp


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: mat_ge_det
    !> @brief Compute determinant of general matrix.
    !> @details
    !! Uses LU factorization to compute the determinant. The input matrix is not changed.
    !----------------------------------------------------------------------------------------------
    function mat_ge_det(a) result(adet)
        real(dp), intent(in) :: a(:, :) !< Matrix.
        real(dp) :: adet !< Determinant.
        real(dp), allocatable :: tmp(:, :) !< Work array.
        integer, allocatable :: ipiv(:)   ! pivot indices
        integer :: sgn !< Sign change based on pivots.
        integer :: info !< LAPACK exit status.
        integer :: i, n

        n = size(a, 1)
        allocate(tmp, source=a)
        allocate(ipiv(n))

        ! LU factorization using partial pivoting with row interchanges.
        call getrf(tmp, ipiv, info)
        if (info /= 0) then
            adet = 0.0_dp
            return
        end if

        sgn = 1
        adet = 1.0_dp
        do i = 1, n
            if (ipiv(i) /= i) sgn = -sgn
            adet = adet * tmp(i, i)
        end do
        adet = sgn * adet
    end function mat_ge_det


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_ge_mmm
    !> @brief Multiply 3 matrices
    !----------------------------------------------------------------------------------------------
    subroutine mat_ge_mmm(a, b, c, d, transa, transc)
        real(dp) :: a(:, :) !< Input matrix 1.
        real(dp) :: b(:, :) !< Input matrix 2.
        real(dp) :: c(:, :) !< Input matrix 3.
        real(dp) :: d(:, :) !< Output matrix.
        character(len=1), intent(in), optional :: transa !< op(A)
        character(len=1), intent(in), optional :: transc !< op(C)

        real(dp), allocatable :: wrk(:, :)
        character(len=1) :: wrk_transa
        character(len=1) :: wrk_transc

        wrk_transa = 'N'
        if (present(transa)) wrk_transa = transa
        wrk_transc = 'N'
        if (present(transc)) wrk_transc = transc

        if (wrk_transa == 'N') then
            allocate(wrk(size(a, 1), size(b, 2)))
        else
            allocate(wrk(size(a, 2), size(b, 2)))
        end if

        call gemm(a, b, wrk, transa=wrk_transa)
        call gemm(wrk, c, d, transb=wrk_transc)
    end subroutine mat_ge_mmm


end module matrix_mod
