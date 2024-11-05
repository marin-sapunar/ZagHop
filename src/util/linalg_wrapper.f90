module linalg_wrapper_mod
    use global_defs
    implicit none
    private

#if BLA_INT64
    integer, parameter :: blas_int = int64
#else
    integer, parameter :: blas_int = int32
#endif

    public :: blas_int
    public :: dot
    public :: gemv
    public :: gemm
    public :: ger
    public :: gesvd
    public :: getrf
    public :: getri
    public :: syev

    interface gemm
        module procedure gemm_d, gemm_z
    end interface gemm

contains


    function dot(x, y) result(z)
        real(dp) :: x(:)
        real(dp) :: y(:)
        real(dp) :: z
        z = sum(x*y)
    end function dot


    subroutine gemv(a, x, y, trans, alpha, beta)
        real(dp) :: a(:, :)
        real(dp) :: x(:)
        real(dp) :: y(:)
        character(len=1), optional :: trans
        real(dp), optional :: alpha
        real(dp), optional :: beta
        character(len=1) :: wrk_trans
        real(dp) :: wrk_alpha
        real(dp) :: wrk_beta
        integer(blas_int) :: m, n
        external dgemv
        m = size(a, 1)
        n = size(a, 2)
        wrk_trans = 'N'
        if (present(trans)) wrk_trans = trans
        wrk_alpha = num1
        if (present(alpha)) wrk_alpha = alpha
        wrk_beta = num0
        if (present(beta)) wrk_beta = beta
        call dgemv(wrk_trans, m, n, wrk_alpha, a, m, x, 1, wrk_beta, y, 1)
    end subroutine gemv


    subroutine gemm_d(a, b, c, transa, transb, alpha, beta)
        real(dp) :: a(:, :)
        real(dp) :: b(:, :)
        real(dp) :: c(:, :)
        character(len=1), optional :: transa
        character(len=1), optional :: transb
        real(dp), optional :: alpha
        real(dp), optional :: beta
        character(len=1) :: wrk_transa
        character(len=1) :: wrk_transb
        real(dp) :: wrk_alpha
        real(dp) :: wrk_beta
        integer(blas_int) :: m, n, k, lda, ldb, ldc
        external dgemm
        lda = size(a, 1)
        ldb = size(b, 1)
        ldc = size(c, 1)
        n = size(c, 2)
        wrk_transa = 'N'
        if (present(transa)) wrk_transa = transa
        wrk_transb = 'N'
        if (present(transb)) wrk_transb = transb
        wrk_alpha = num1
        if (present(alpha)) wrk_alpha = alpha
        wrk_beta = num0
        if (present(beta)) wrk_beta = beta
        m = lda
        if (wrk_transa /= 'N') then
            m = size(a, 2)
            k = lda
        else
            m = lda
            k = size(a, 2)
        end if
        call dgemm(wrk_transa, wrk_transb, m, n, k, wrk_alpha, a, lda, b, ldb, wrk_beta, c, ldc)
    end subroutine gemm_d


    subroutine gemm_z(a, b, c, transa, transb, alpha, beta)
        complex(dp) :: a(:, :)
        complex(dp) :: b(:, :)
        complex(dp) :: c(:, :)
        character(len=1), optional :: transa
        character(len=1), optional :: transb
        complex(dp), optional :: alpha
        complex(dp), optional :: beta
        character(len=1) :: wrk_transa
        character(len=1) :: wrk_transb
        complex(dp) :: wrk_alpha
        complex(dp) :: wrk_beta
        integer(blas_int) :: m, n, k, lda, ldb, ldc
        external zgemm
        lda = size(a, 1)
        ldb = size(b, 1)
        ldc = size(c, 1)
        n = size(c, 2)
        wrk_transa = 'N'
        if (present(transa)) wrk_transa = transa
        wrk_transb = 'N'
        if (present(transb)) wrk_transb = transb
        wrk_alpha = cmplx(num1, num0, dp)
        if (present(alpha)) wrk_alpha = alpha
        wrk_beta = cmplx(num0, num0, dp)
        if (present(beta)) wrk_beta = beta
        m = lda
        if (wrk_transa /= 'N') then
            m = size(a, 2)
            k = lda
        else
            m = lda
            k = size(a, 2)
        end if
        call zgemm(wrk_transa, wrk_transb, m, n, k, wrk_alpha, a, lda, b, ldb, wrk_beta, c, ldc)
    end subroutine gemm_z


    subroutine ger(a, x, y, alpha)
        real(dp) :: x(:)
        real(dp) :: y(:)
        real(dp) :: a(:, :)
        real(dp), optional :: alpha
        integer(blas_int) :: m
        integer(blas_int) :: n
        real(dp) :: wrk_alpha
        external dger
        m = size(x)
        n = size(y)
        wrk_alpha = num1
        if (present(alpha)) wrk_alpha = alpha
        call dger(m, n, wrk_alpha, x, 1, y, 1, a, m)
    end subroutine ger


    subroutine gesvd(a, s, u, vt, job, info)
        real(dp) :: a(:, :)
        real(dp) :: s(:)
        real(dp), optional :: u(:, :)
        real(dp), optional :: vt(:, :)
        character(len=1), optional :: job
        integer, optional :: info
        real(dp), allocatable :: wrk_u(:, :)
        real(dp), allocatable :: wrk_vt(:, :)
        real(dp), allocatable :: work(:)
        character(len=1) :: wrk_job
        character(len=1) :: jobu
        character(len=1) :: jobvt
        integer(blas_int) :: wrk_info
        integer(blas_int) :: m, n, ldvt, ldu, lwork
        external :: dgesvd

        m = size(a, 1)
        n = size(a, 2)
        wrk_job = 'N'
        if (present(job)) wrk_job = job
        if (present(u)) then
            allocate(wrk_u, source=u)
            ldu = size(wrk_u, 1)
            if (size(u, 2) == m) then
                jobu = 'A'
            else
                jobu = 'S'
            end if
        else
            if (wrk_job == 'U') then
                jobu = 'O'
            else
                jobu = 'N'
            end if
        end if
        if (present(vt)) then
            allocate(wrk_vt, source=vt)
            ldvt = size(wrk_vt, 1)
            if (size(vt, 1) == n) then
                jobvt = 'A'
            else
                jobvt = 'S'
            end if
        else
            if (wrk_job == 'V') then
                jobvt = 'O'
            else
                jobvt = 'N'
            end if
        end if

        ! Determine size of work array.
        allocate(work(1))
        lwork = -1
        call dgesvd(jobu, jobvt, m, n, a, m, s, wrk_u, ldu, wrk_vt, ldvt, work, lwork, wrk_info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Call DGESVD.
        call dgesvd(jobu, jobvt, m, n, a, m, s, wrk_u, ldu, wrk_vt, ldvt, work, lwork, wrk_info)
        if (present(u)) u = wrk_u
        if (present(vt)) vt = wrk_vt
        if (present(info)) then
            info = wrk_info
        else
            if (wrk_info /= 0) then
                write(stderr, *) 'Error. GESVD call failed.'
                write(stderr, '(a, i0)') '        Info: ', wrk_info
                stop
            end if
        end if
    end subroutine gesvd


    subroutine getrf(a, ipiv, info)
        real(dp) :: a(:, :)
        integer, optional :: ipiv(:)
        integer, optional :: info
        integer(blas_int) :: wrk_info
        integer(blas_int), allocatable :: wrk_ipiv(:)
        integer(blas_int) :: n, m
        external :: dgetrf

        m = size(a, 1)
        n = size(a, 2)
        if (present(ipiv)) then
            allocate(wrk_ipiv(size(ipiv)))
        else
            allocate(wrk_ipiv(min(m, n)))
        end if

        ! Call DGETRF
        call dgetrf(m, n, a, n, wrk_ipiv, wrk_info)
        if (present(ipiv)) ipiv = wrk_ipiv
        if (present(info)) then
            info = wrk_info
        else
            if (wrk_info /= 0) then
                write(stderr, *) 'Error. GETRF call failed.'
                write(stderr, '(a, i0)') '        Info: ', wrk_info
                stop
            end if
        end if
    end subroutine getrf

    subroutine getri(a, ipiv, info)
        real(dp) :: a(:, :)
        integer  :: ipiv(:)
        integer, optional :: info
        integer(blas_int) :: lwrk
        integer(blas_int), allocatable :: wrk(:)
        integer(blas_int) :: n, lda
        external :: dgetri

        lda = size(a, 1)
        n = size(a, 2)
        lwrk = n
        allocate(wrk(lwrk))

        ! Call DGETRI
        call  dgetri( n, a, lda, ipiv, wrk, lwrk, info )
        if (info /= 0) then
            write(stderr, *) 'Error. GETRI call failed.'
            write(stderr, '(a, i0)') '        Info: ', info
            stop
        end if
    end subroutine getri

    subroutine syev(a, w, jobz, uplo, info)
        real(dp) :: a(:, :)
        real(dp) :: w(:)
        character(len=1), optional :: jobz
        character(len=1), optional :: uplo
        integer, optional :: info
        real(dp), allocatable :: work(:)
        character(len=1) :: wrk_jobz
        character(len=1) :: wrk_uplo
        integer(blas_int) :: wrk_info
        integer(blas_int) :: n, lwork
        external :: dsyev

        n = size(a, 1)
        wrk_jobz = 'N'
        if (present(jobz)) wrk_jobz = jobz
        wrk_uplo = 'U'
        if (present(uplo)) wrk_uplo = uplo

        ! Determine size of work array.
        allocate(work(1))
        lwork = -1
        call dsyev(wrk_jobz, wrk_uplo, n, a, n, w, work, lwork, wrk_info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Call DSYEV.
        call dsyev(wrk_jobz, wrk_uplo, n, a, n, w, work, lwork, wrk_info)
        if (present(info)) then
            info = wrk_info
        else
            if (wrk_info /= 0) then
                write(stderr, *) 'Error. SYEV call failed.'
                write(stderr, '(a, i0)') '        Info: ', wrk_info
                stop
            end if
        end if
    end subroutine syev


end module linalg_wrapper_mod
