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
    type, extends(rng_type_int) :: rng_pcg_xsh_rr
        logical :: initialized = .false.
        integer(int64) :: state
    contains
        procedure :: init => xsh_rr_init
        procedure :: random_integer => xsh_rr_int
        procedure :: save_state => xsh_rr_save_state
        procedure :: load_state => xsh_rr_load_state
    end type rng_pcg_xsh_rr


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_init
    !> @brief Initialize the random number generator.
    !> @details
    !! Calls the init of the parent class, and initializes the state vector based on the seed.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_init(self, seed)
        class(rng_pcg_xsh_rr), intent(inout) :: self
        integer, intent(in) :: seed !< Seed used for the initialization. 

        call self%init_seed(seed)
        call set_range(self)
        self%state = self%seed + incr
        self%state = self%state * mult + incr !< Advance state once.
        self%initialized = .true.
    end subroutine xsh_rr_init


    !> Range of the integers returned by the generator, for converting them to reals.
    subroutine set_range(self)
        class(rng_pcg_xsh_rr), intent(inout) :: self

        self%irange = 2.0**32_dp
        self%i_irangem1 = 1.0_dp / (self%irange - 1.0_dp)
    end subroutine set_range


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_save_state
    !> @brief Add the state of the generator to a JSON object, in the object pcg_xsh_rr.
    !> @details
    !! The 64-bit state is written as a decimal string, since JSON integers of json-fortran are
    !! 32-bit.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(rng_pcg_xsh_rr), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: pcg
        character(len=24) :: state_str

        call rng_save_state(self, json, obj)
        call json%create_object(pcg, 'pcg_xsh_rr')
        call json%add(pcg, 'initialized', self%initialized)
        write(state_str, '(i0)') self%state
        call json%add(pcg, 'state', trim(state_str))
        call json%add(obj, pcg)
    end subroutine xsh_rr_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_load_state
    !> @brief Restore the state of the generator saved by save_state.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child
        class(rng_pcg_xsh_rr), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: pcg
        character(len=:), allocatable :: state_str
        integer :: ios

        call rng_load_state(self, json, obj)
        pcg => json_get_child(json, obj, 'pcg_xsh_rr')
        call json_get_required(json, pcg, 'initialized', self%initialized)
        call json_get_required(json, pcg, 'state', state_str)
        read(state_str, *, iostat=ios) self%state
        if (ios /= 0) then
            write(stderr, *) 'Error in random_pcg_mod, xsh_rr_load_state subroutine.'
            write(stderr, *) '  Cannot read the state of the generator: ', state_str
            stop 1
        end if
        call set_range(self)
    end subroutine xsh_rr_load_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xsh_rr_int
    !> @brief Get a random integer from the generator.
    !----------------------------------------------------------------------------------------------
    subroutine xsh_rr_int(self, rnum)
        class(rng_pcg_xsh_rr), intent(inout) :: self
        integer(int32), intent(out) :: rnum
        integer(int64) :: x
        integer(int32) :: count

        ! Initialize the generator if it hasn't been done before.
        if (.not. self%initialized) call self%init(-1)

        x = self%state
        count = int(shiftr(x, 59), kind=int32)
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
