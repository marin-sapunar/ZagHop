
#if __INTEL_COMPILER
include 'mkl_vsl.f90'
#endif

module random_mod
    use global_defs
#if __INTEL_COMPILER
    use mkl_vsl_type
    use mkl_vsl
#endif
    implicit none

    private
    public :: init_random_seed
    public :: random_gaussian

#if __INTEL_COMPILER
    type(VSL_STREAM_STATE) :: rng_stream
#endif


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init_random_seed
    !> @brief Initialize the random number generator.
    !> @details
    !! If a positive integer is passed to the subroutine it is used as the seed. Otherwise, a seed
    !! is generated automatically using the current time and PID of the process. This routine
    !! initializes both the native Fortran random_number subroutine and an Intel Vector Statistics
    !! stream rng_stream.
    !
    !> @note This module should be compiled using a C-style preprocessor. For the intel compiler,
    !! this is set using the -fpp flag.
    !----------------------------------------------------------------------------------------------
    subroutine init_random_seed(tseed)
#if __INTEL_COMPILER
use ifport
#endif
        integer, intent(inout) :: tseed !< Seed used for the initialization. 
        integer :: tclock
        integer, allocatable :: aseed(:)
        integer :: n
        integer :: i
        integer :: pid
        real(dp) :: rnum(50)

        call random_seed(size = n)
        allocate(aseed(n))
        if (tseed >= 0) then
            aseed = tseed + 1307 * [(i - 1, i = 1, n)]
            call random_seed(put = aseed)
        else
            call system_clock(count = tclock)
            pid = getpid()
            aseed = tclock + 271 * pid + 1307 * [(i - 1, i = 1, n)]
            call random_seed(put = aseed)
            tseed = aseed(1)
        end if
        ! When running multiple programs in the background, the time is usually identical and the
        ! pid is different by 1. In this case, (for ifort compiled programs) the first few random
        ! numbers are almost identical so we move the generator further along.
        call random_number(rnum)

#if __INTEL_COMPILER
        n = vslnewstream(rng_stream, VSL_BRNG_MCG31, tseed)
        if (n /= VSL_STATUS_OK) then
            write(stderr, *) "Error in random_mod, init_random_seed subroutine."
            write(stderr, *) " vslnewstream exit code:", n
            call abort()
        end if
#endif
    end subroutine init_random_seed


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: random_gaussian
    !> @brief Sample random numbers from a Gaussian distribution.
    !> @details
    !! Wrapper for the vdrnggaussian subroutine from MKL. Assumes that the init_random_seed
    !! subroutine has been called earlier.
    !----------------------------------------------------------------------------------------------
    subroutine random_gaussian(mean, sigma, rnum)
        real(dp), intent(in) :: mean
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: rnum(:)
        complex(dp)           :: gpair
        integer :: method
        integer :: n
        integer(i4) :: err

#if __INTEL_COMPILER
        method = VSL_RNG_METHOD_GAUSSIAN_ICDF
        n = size(rnum)
        err = vdrnggaussian(method, rng_stream, n, rnum, mean, sigma)

        if (err /= VSL_STATUS_OK) then
            write(stderr, *) "Error in random_mod, rnggaussian subroutine."
            write(stderr, *) " vdrnggaussian exit code:", n
            call abort()
        end if
#elseif QUANTICS
       do i=1,n
          gpair=randgaussian()
          rnum(i)=mean + sigma*dble(gpair)
       enddo
#endif
    end subroutine random_gaussian


end module random_mod
