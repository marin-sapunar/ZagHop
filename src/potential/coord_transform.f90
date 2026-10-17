!--------------------------------------------------------------------------------------------------
! MODULE: coord_transform_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2026
!
! DESCRIPTION:
!> @brief Transformations between the coordinate systems used for the classical propagation and
!! for the evaluation of the potential.
!--------------------------------------------------------------------------------------------------
module coord_transform_mod
    use global_defs
    use file_mod, only : reader
    implicit none

    private
    public :: normal_modes

    !> Dalton to electron rest mass. Same as in the constants module, which is not part of the
    !! potential library.
    real(dp), parameter :: Da_me = 1822.88849_dp


    !----------------------------------------------------------------------------------------------
    ! TYPE: nmode
    !> @brief Normal modes at a reference geometry.
    !> @details
    !! Dimensionless normal mode coordinates are Q = Omega^(1/2) L^T M^(1/2) (R - R0), where L are
    !! the mass-weighted normal modes, Omega the frequencies, M the masses and R0 the reference
    !! geometry. Modes with zero frequency are excluded from Q.
    !----------------------------------------------------------------------------------------------
    type normal_modes
        character(len=:), allocatable :: v0_file
        integer :: natom = 0
        integer :: nfreq = 0 !< Number of frequencies, including zero frequencies.
        integer :: zero_mode = 0 !< Number of zero frequencies.
        character(len=2), allocatable :: sym(:) !< Atom symbols.
        integer, allocatable :: atnum(:) !< Atomic numbers.
        real(dp), allocatable :: mass(:) !< Atom masses (amu).
        real(dp), allocatable :: ref_geom(:, :) !< Reference geometry (bohr), shape (3, natom).
        real(dp), allocatable :: freq(:) !< Frequencies, shape (nfreq).
        real(dp), allocatable :: nm(:, :) !< Mass-weighted normal modes, shape (3*natom, nfreq).
        real(dp), allocatable :: c_to_q(:, :) !< Cartesian displacement to Q transformation,
                                              !! shape (nfreq-zero_mode, 3*natom).
    contains
        procedure :: read_v0
        procedure :: cart_to_nm
        procedure :: deriv_nm_to_cart
    end type normal_modes


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_v0
    ! DESCRIPTION:
    !> @brief Read reference geometry, frequencies and normal modes from a V0 file.
    !> @details
    !! The geometry and normal modes are read from the "Geometry" and "Mass-weighted normal modes"
    !! sections. If either section is missing, a warning is printed and the transformation is not
    !! available.
    !----------------------------------------------------------------------------------------------
    subroutine read_v0(self, v0_file)
        class(normal_modes) :: self
        character(len=*), intent(in) :: v0_file
        type(reader) :: readf
        logical :: have_line
        integer :: i, ncol, nrow
        character(len=2) :: at_sym
        real(dp) :: at_num, at_mass, xyz(3)
        real(dp), allocatable :: row(:), wrk(:), sqrt_mass(:)

        self%v0_file = v0_file
        call readf%open(v0_file, abort_on_eof=.false.)

        ! Geometry and normal mode sections have no explicit length. They end at the first line
        ! that does not match the section format, and that line is processed in the next cycle.
        have_line = .false.
        do
            if (.not. have_line) call readf%next()
            have_line = .false.
            if (is_iostat_end(readf%iostat)) exit
            call readf%parseline(' ')
            select case(readf%args(1)%s)
            case('Geometry')
                self%sym = [character(len=2) ::]
                self%atnum = [integer ::]
                self%mass = [real(dp) ::]
                wrk = [real(dp) ::]
                do
                    call readf%next()
                    if (is_iostat_end(readf%iostat)) exit
                    call readf%parseline(' ')
                    if (readf%narg /= 6) exit
                    read(readf%line, *) at_sym, at_num, xyz, at_mass
                    self%sym = [self%sym, at_sym]
                    self%atnum = [self%atnum, nint(at_num)]
                    self%mass = [self%mass, at_mass]
                    wrk = [wrk, xyz]
                end do
                self%natom = size(self%sym)
                self%ref_geom = reshape(wrk, [3, self%natom])
                have_line = .true.
            case('Frequencies')
                call readf%next()
                call readf%parseline(' ')
                self%nfreq = readf%narg
                allocate(self%freq(self%nfreq))
                read(readf%line, *) self%freq
            case('Mass-weighted')
                wrk = [real(dp) ::]
                ncol = 0
                nrow = 0
                do
                    call readf%next()
                    if (is_iostat_end(readf%iostat)) exit
                    call readf%parseline(' ')
                    if (nrow == 0) then
                        ncol = readf%narg
                        allocate(row(ncol))
                    end if
                    if (readf%narg /= ncol) exit
                    read(readf%line, *) row
                    wrk = [wrk, row]
                    nrow = nrow + 1
                end do
                self%nm = transpose(reshape(wrk, [ncol, nrow]))
                have_line = .true.
            end select
        end do
        call readf%close()

        if (.not. allocated(self%ref_geom)) then
            write(stderr, *) 'Warning in coord_transform_mod, read_v0 subroutine.'
            write(stderr, *) '  Geometry section not found in file: ', v0_file
        end if
        if (.not. allocated(self%nm)) then
            write(stderr, *) 'Warning in coord_transform_mod, read_v0 subroutine.'
            write(stderr, *) '  Mass-weighted normal modes section not found in file: ', v0_file
        end if
        if (.not. (allocated(self%ref_geom) .and. allocated(self%nm))) return

        if (size(self%nm, 1) /= 3 * self%natom) then
            write(stderr, *) 'Error in coord_transform_mod, read_v0 subroutine.'
            write(stderr, *) '  Number of rows in normal mode matrix does not match 3*natom.'
            stop
        end if
        if (size(self%nm, 2) /= self%nfreq) then
            write(stderr, *) 'Error in coord_transform_mod, read_v0 subroutine.'
            write(stderr, *) '  Number of columns in normal mode matrix does not match'
            write(stderr, *) '  number of frequencies.'
            stop
        end if

        self%zero_mode = count(abs(self%freq) < 1.0e-12_dp)
        sqrt_mass = reshape(spread(sqrt(self%mass * Da_me), 1, 3), [3 * self%natom])
        allocate(self%c_to_q(self%nfreq - self%zero_mode, 3 * self%natom))
        do i = 1, self%nfreq - self%zero_mode
            self%c_to_q(i, :) = sqrt(self%freq(self%zero_mode + i))                              &
            &                   * self%nm(:, self%zero_mode + i) * sqrt_mass
        end do
    end subroutine read_v0


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: cart_to_nm
    ! DESCRIPTION:
    !> @brief Convert Cartesian geometry to dimensionless normal mode coordinates.
    !----------------------------------------------------------------------------------------------
    function cart_to_nm(self, geom) result(q)
        class(normal_modes), intent(in) :: self
        real(dp), intent(in) :: geom(:, :) !< Cartesian geometry (bohr), shape (3, natom).
        real(dp), allocatable :: q(:)

        call check_available(self, 'cart_to_nm')
        if (size(geom, 1) /= 3 .or. size(geom, 2) /= self%natom) then
            write(stderr, *) 'Error in coord_transform_mod, cart_to_nm function.'
            write(stderr, '(3x,a,i0,a)') 'Geometry should have dimensions (3, ', self%natom, ').'
            stop
        end if

        q = matmul(self%c_to_q, reshape(geom - self%ref_geom, [3 * self%natom]))
    end function cart_to_nm


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: deriv_nm_to_cart
    ! DESCRIPTION:
    !> @brief Convert derivatives with respect to normal modes to Cartesian derivatives.
    !> @details
    !! Applies to any derivative with respect to the nuclear coordinates (gradients, nonadiabatic
    !! coupling vectors): dX/dR = M^(1/2) L Omega^(1/2) dX/dQ.
    !----------------------------------------------------------------------------------------------
    function deriv_nm_to_cart(self, dq) result(dr)
        class(normal_modes), intent(in) :: self
        real(dp), intent(in) :: dq(:) !< Derivative with respect to Q.
        real(dp), allocatable :: dr(:) !< Derivative with respect to R, shape (3*natom).

        call check_available(self, 'deriv_nm_to_cart')
        if (size(dq) /= size(self%c_to_q, 1)) then
            write(stderr, *) 'Error in coord_transform_mod, deriv_nm_to_cart function.'
            write(stderr, '(3x,a,i0,a)') 'Derivative should have dimension ',                  &
            &                            size(self%c_to_q, 1), '.'
            stop
        end if

        dr = matmul(dq, self%c_to_q)
    end function deriv_nm_to_cart


    subroutine check_available(self, routine)
        class(normal_modes), intent(in) :: self
        character(len=*), intent(in) :: routine

        if (allocated(self%c_to_q)) return
        write(stderr, *) 'Error in coord_transform_mod, '//routine//' function.'
        write(stderr, *) '  Transformation between Cartesian and normal mode coordinates requires'
        write(stderr, *) '  the reference geometry and normal modes from the V0 file.'
        if (allocated(self%v0_file)) write(stderr, *) '  File: ', self%v0_file
        stop
    end subroutine check_available


end module coord_transform_mod
