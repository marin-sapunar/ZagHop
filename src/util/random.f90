!--------------------------------------------------------------------------------------------------
! MODULE: random_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2024
!
! DESCRIPTION:
!> @brief Module containing the `rng_type` base class.
!--------------------------------------------------------------------------------------------------
module random_mod
    use global_defs
    implicit none

    private
    public :: rng_type
    public :: rng_type_int


    !----------------------------------------------------------------------------------------------
    ! TYPE: rng_type
    !> @brief Random number generator base class.
    !> @details
    !! This is the base class which defines helper procedures to be used along with the actual RNG
    !! which are implemented as child classes.
    !----------------------------------------------------------------------------------------------
    type, abstract :: rng_type
        integer :: seed = -1 !< Seed used for initializing the RNG.
    contains
        private
        procedure, public :: init_seed => rng_init_seed
        procedure(rngReal), deferred, public :: random_real
        procedure(rngInt), deferred, public :: random_integer
        procedure(rngInit), deferred, public :: init
        procedure :: uniform_1d => rng_uniform_1d
        procedure :: uniform_2d => rng_uniform_2d
        procedure :: uniform_3d => rng_uniform_3d
        procedure :: integer_1d => rng_integer_1d
        procedure :: integer_2d => rng_integer_2d
        procedure :: integer_3d => rng_integer_3d
        procedure :: random_gaussian => rng_gaussian
        procedure :: gaussian_1d => rng_gaussian_1d
        procedure :: gaussian_2d => rng_gaussian_2d
        procedure :: gaussian_3d => rng_gaussian_3d
        generic, public :: uniform => random_real, uniform_1d, uniform_2d, uniform_3d
        generic, public :: integer => random_integer, integer_1d, integer_2d, integer_3d
        generic, public :: gaussian => random_gaussian, gaussian_1d, gaussian_2d, gaussian_3d
    end type


    !----------------------------------------------------------------------------------------------
    ! TYPE: rng_type_int
    !> @brief Base class for RNGs which implement the random_integer method.
    !> @details
    !! This is the base class which assumes that child classes will implement the `random_integer`
    !! method along with the `init` method. The `random_real` method here will convert the outputs
    !! of `random_integer` to real random numbers between 0 and 1.
    !----------------------------------------------------------------------------------------------
    type, abstract, extends(rng_type) :: rng_type_int
        real(dp) :: irange = 0.0_dp !< Range of integers returned by the RNG.
        real(dp) :: i_irangem1 = 0.0_dp !< 1 / (irange - 1)
    contains
        procedure :: random_real => rng_real
    end type


    abstract interface
        subroutine rngReal(self, rnum)
            import rng_type
            import dp
            class(rng_type), intent(inout) :: self
            real(dp), intent(out) :: rnum
        end subroutine rngReal
    end interface


    abstract interface
        subroutine rngInt(self, rnum)
            import rng_type
            import int32
            class(rng_type), intent(inout) :: self
            integer(int32), intent(out) :: rnum
        end subroutine rngInt
    end interface


    abstract interface
        subroutine rngInit(self, seed)
            import rng_type
            class(rng_type), intent(inout) :: self
            integer, intent(in) :: seed
        end subroutine rngInit
    end interface


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_init_seed
    !> @brief Initialize a seed for the random number generator.
    !> @details
    !! If a positive integer is passed to the subroutine it is used as the seed. Otherwise, a seed
    !! is generated automatically using the current time and PID of the process.
    !----------------------------------------------------------------------------------------------
    subroutine rng_init_seed(self, seed)
#if __INTEL_COMPILER
        use ifport, only : getpid
#endif
        class(rng_type) :: self
        integer, intent(in) :: seed !< Seed used for the initialization. 
        integer :: tclock
        integer :: pid

        if (seed <= 0) then
            call system_clock(count = tclock)
            pid = getpid()
            self%seed = tclock + 271 * pid
        else
            self%seed = seed
        end if
    end subroutine rng_init_seed


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_real
    !> @brief Return a random real number between 0 and 1.
    !----------------------------------------------------------------------------------------------
    subroutine rng_real(self, rnum)
        class(rng_type_int), intent(inout) :: self
        real(dp), intent(out) :: rnum
        integer(int32) :: z

        call self%random_integer(z)
        if (z < 0) then
            rnum = (real(z, kind=dp) + self%irange) * self%i_irangem1
        else
            rnum = real(z, kind=dp) * self%i_irangem1
        end if
    end subroutine rng_real


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_uniform_1d
    !> @brief Fill a 1D array with uniform random real numbers between 0 and 1.
    !----------------------------------------------------------------------------------------------
    subroutine rng_uniform_1d(self, rnum)
        class(rng_type) :: self
        real(dp), intent(out) :: rnum(:)
        integer :: i

        do i = 1, size(rnum)
            call self%uniform(rnum(i))
        end do
    end subroutine rng_uniform_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_uniform_2d
    !> @brief Fill a 2D array with uniform random real numbers between 0 and 1.
    !----------------------------------------------------------------------------------------------
    subroutine rng_uniform_2d(self, rnum)
        class(rng_type) :: self
        real(dp), intent(out) :: rnum(:, :)
        integer :: i

        do i = 1, size(rnum, 2)
            call self%uniform(rnum(:, i))
        end do
    end subroutine rng_uniform_2d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_uniform_3d
    !> @brief Fill a 3D array with uniform random real numbers between 0 and 1.
    !----------------------------------------------------------------------------------------------
    subroutine rng_uniform_3d(self, rnum)
        class(rng_type) :: self
        real(dp), intent(out) :: rnum(:, :, :)
        integer :: i

        do i = 1, size(rnum, 3)
            call self%uniform(rnum(:, :, i))
        end do
    end subroutine rng_uniform_3d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_integer_1d
    !> @brief Fill a 1D array with random integers.
    !----------------------------------------------------------------------------------------------
    subroutine rng_integer_1d(self, rnum)
        class(rng_type) :: self
        integer(int32), intent(out) :: rnum(:)
        integer :: i

        do i = 1, size(rnum)
            call self%integer(rnum(i))
        end do
    end subroutine rng_integer_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_integer_2d
    !> @brief Fill a 2D array with random integers.
    !----------------------------------------------------------------------------------------------
    subroutine rng_integer_2d(self, rnum)
        class(rng_type) :: self
        integer(int32), intent(out) :: rnum(:, :)
        integer :: i

        do i = 1, size(rnum, 2)
            call self%integer(rnum(:, i))
        end do
    end subroutine rng_integer_2d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_integer_3d
    !> @brief Fill a 3D array with random integers.
    !----------------------------------------------------------------------------------------------
    subroutine rng_integer_3d(self, rnum)
        class(rng_type) :: self
        integer(int32), intent(out) :: rnum(:, :, :)
        integer :: i

        do i = 1, size(rnum, 3)
            call self%integer(rnum(:, :, i))
        end do
    end subroutine rng_integer_3d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: random_gaussian
    !> @brief Sample random numbers from a Gaussian distribution.
    !> @details
    !! Uses a polar Box-Muller transform to generate pairs of independent normally distributed
    !! random numbers.
    !----------------------------------------------------------------------------------------------
    subroutine rng_gaussian(self, mean, sigma, rnum)
        class(rng_type) :: self
        real(dp), intent(in) :: mean
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: rnum
        logical, save :: flag = .false.
        real(dp), save :: gpair(2) = 0.0_dp
        real(dp) :: rsq

        if (flag) then
            rnum = gpair(2)
            flag = .false.
        else
            do
                call self%uniform(gpair)
                gpair = 2.0_dp * gpair - 1.0_dp
                rsq = sum(gpair**2)
                if ((rsq <= 1.0_dp) .and. (rsq /= 0.0_dp)) exit
            end do
            gpair = sqrt(-2.0_dp * log(rsq) / rsq) * gpair
            rnum = gpair(1)
            flag = .true.
        end if

        rnum = mean + sigma * rnum
    end subroutine rng_gaussian


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_gaussian_1d
    !> @brief Sample random numbers from a Gaussian distribution.
    !----------------------------------------------------------------------------------------------
    subroutine rng_gaussian_1d(self, mean, sigma, rnum)
        class(rng_type) :: self
        real(dp), intent(in) :: mean
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: rnum(:)
        integer :: i

        do i = 1, size(rnum)
            call self%gaussian(mean, sigma, rnum(i))
        end do
    end subroutine rng_gaussian_1d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_gaussian_2d
    !> @brief Sample random numbers from a Gaussian distribution.
    !----------------------------------------------------------------------------------------------
    subroutine rng_gaussian_2d(self, mean, sigma, rnum)
        class(rng_type) :: self
        real(dp), intent(in) :: mean
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: rnum(:, :)
        integer :: i

        do i = 1, size(rnum, 2)
            call self%gaussian(mean, sigma, rnum(:, i))
        end do
    end subroutine rng_gaussian_2d


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rng_gaussian_3d
    !> @brief Sample random numbers from a Gaussian distribution.
    !----------------------------------------------------------------------------------------------
    subroutine rng_gaussian_3d(self, mean, sigma, rnum)
        class(rng_type) :: self
        real(dp), intent(in) :: mean
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: rnum(:, :, :)
        integer :: i

        do i = 1, size(rnum, 3)
            call self%gaussian(mean, sigma, rnum(:, :, i))
        end do
    end subroutine rng_gaussian_3d


end module random_mod
