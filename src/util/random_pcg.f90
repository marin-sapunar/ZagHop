!--------------------------------------------------------------------------------------------------
! MODULE: random_pcg_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2024
!
! DESCRIPTIONi:
!> @brief Module containing the `rng_pcg_xsh_rr` type.
!> @details
!! Implementation of the Permuted congruential generator variant PCG-XSH-RR by M. O'Neill with a
!! period of 2^64 and 32-bit output. The PCG RNGs are described in the technical report "PCG: A 
!! Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number
!! Generation" found at https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf. The RNG code here is a
!! translation of the C++ code from the PCG Wikipedia page.
!--------------------------------------------------------------------------------------------------
module random_pcg_mod
    use global_defs
    use random_mod
    implicit none


    private
    public :: rng_pcg_xsh_rr


    ! Parameters of the RNG.
    integer(int64), parameter :: mult = int(Z'5851F42D4C957F2D', kind=int64)
    integer(int64), parameter :: incr = int(Z'14057B7EF767814F', kind=int64)
    

    !----------------------------------------------------------------------------------------------
    ! TYPE: rng_pcg_xsh_rr
    !> @brief Permuted congruential generator variant PCG-XSH-RR.
    !----------------------------------------------------------------------------------------------
    type, extends(rng_type) :: rng_pcg_xsh_rr
        logical :: initialized = .false.
        integer(int64) :: state
    contains
        procedure :: init => xsh_rr_init
        procedure :: random_int => xsh_rr_int
    end type rng_pcg_xsh_rr


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_init
    !> @brief Initialize the random number generator.
    !> @details
    !! Calls the init of the parent class, and initializes the state vector based on the seed.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_init(self, seed)
        class(rng_pcg_xsh_rr) :: self
        integer, intent(in) :: seed !< Seed used for the initialization. 
        integer :: i
        integer(int32) :: z

        call self%rng_type%init(seed)
        self%rng_out = 1
        self%irange = 2.0**32_dp
        self%i_irangem1 = 1.0_dp / (self%irange - 1.0_dp)
        self%state = self%seed + incr
        self%state = self%state * mult + incr !< Advance state once.
        self%initialized = .true.
    end subroutine xsh_rr_init


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_int
    !> @brief Get a random integer from the generator.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_int(self, rnum)
        class(rng_pcg_xsh_rr) :: self
        integer(int32), intent(out) :: rnum
        integer(int64) :: x
        integer(int32) :: count

        ! Initialize the generator if it hasn't been done before.
        if (.not. self%initialized) call self%init(-1)

        x = self%state
        count = shiftr(x, 59)
        self%state = x * mult + incr
        x = ieor(x, shiftr(x, 18))
        rnum = rotr32(int((shiftr(x, 27)), kind=int32), count)
    end subroutine xsh_rr_int


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: xsh_rr_int
    !> @brief 
    !----------------------------------------------------------------------------------------------
    function rotr32(x, r) result(rotx)
       integer(int32), intent(in) :: x 
       integer(int32), intent(in) :: r
       integer(int32) :: rotx
       rotx = ior(shiftr(x, r), shiftl(x, iand(-r, 31_int32)))
    end function rotr32


end module random_pcg_mod
