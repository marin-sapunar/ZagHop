module random_mod
    implicit none

    private
    public :: init_random_seed


contains

    
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init_random_seed
    !> @brief Initialize the random number generator.
    !> @details
    !! If a positive integer is passed to the subroutine it is used as the seed. Otherwise, a seed
    !! is generated automatically using the current time and PID of the process.
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
        real :: rnum(50)

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
    end subroutine init_random_seed


end module random_mod
