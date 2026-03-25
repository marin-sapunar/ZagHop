!--------------------------------------------------------------------------------------------------
! MODULE: vc_model_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date January, 2026
!
! DESCRIPTION: 
!> @brief Simple vibronic coupling interface.
!--------------------------------------------------------------------------------------------------
module vc_model_mod
    use global_defs
    use string_mod
    use vc_block_mod
    use file_mod, only : reader
    implicit none

    private
    public :: vc_model

    type vc_model
        character(len=:), allocatable :: template_file
        character(len=:), allocatable :: v0_file
        integer :: nstate(3) = [0, 0, 0]
        integer :: nmode = 0
        integer :: tot_ns = 0
        integer :: zero_mode = 0
        real(dp), allocatable :: freq(:)
        type(vc_block) :: w(3)
        type(vc_block) :: soc
        type(vc_block) :: dm(3) 
        type(vc_block), allocatable :: dw(:, :)
    contains
        procedure :: init => vc_init
        procedure :: read_v0
    end type vc_model


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: vc_init
    !
    ! DESCRIPTION:
    !> @brief Initialize LVC model from template file.
    !> @details
    !----------------------------------------------------------------------------------------------
    subroutine vc_init(self, template_file)
        class(vc_model) :: self
        character(len=*), intent(in) :: template_file
        type(reader) :: readf
        integer :: n_val, i, j, k, mult
        integer ::  ist1, ist2, imode1, imode2
        integer :: slash_pos
        character(len=:), allocatable :: template_dir
        real(dp) :: val
        real(dp), allocatable :: wrk(:, :)
        character(len=3) :: sec
        character(len=1) :: part

        self%template_file = template_file

        slash_pos = index(template_file, '/', back=.true.)
        if (slash_pos > 0) then
            template_dir = template_file(1:slash_pos)
        else
            template_dir = ''
        end if

        call readf%open(self%template_file, abort_on_eof=.false.)

        call readf%next()
        call readf%parseline(' ')
        self%v0_file = template_dir // readf%args(1)%s
        call self%read_v0(self%v0_file)

        call readf%next()
        call readf%parseline(' ')
        if (readf%narg /= 3) then
            write(stderr, *) 'Error in vibronic_mod, vc_init subroutine.'
            write(stderr, *) '  Expected 3 arguments on second line of template file:'
            write(stderr, *) '    n_singlets n_doublets n_triplets'
            stop
        end if
        read(readf%args(1)%s, *) self%nstate(1)
        read(readf%args(2)%s, *) self%nstate(2)
        read(readf%args(3)%s, *) self%nstate(3)

        self%tot_ns = sum(self%nstate * [1, 2, 3])

        do i = 1, 3
            call self%w(i)%allocate(2, self%nmode, self%nstate(i))
            do j = 1, self%nstate(i)
                do k = 1, self%nmode
                    self%w(i)%quad(k, k, j, j) = 0.5_dp * self%freq(k + self%zero_mode)
                end do
            end do
        end do

        allocate(wrk(self%tot_ns, self%tot_ns))

        do
            call readf%next()
            if (is_iostat_end(readf%iostat)) exit
            call readf%parseline(' ')
            select case(readf%args(1)%s)
            case('epsilon')
                call readf%next()
                call readf%parseline(' ')
                read(readf%args(1)%s, *) n_val
                do i = 1, n_val
                    call readf%next()
                    call readf%parseline(' ')
                    read(readf%args(1)%s, *) mult
                    read(readf%args(2)%s, *) ist1
                    read(readf%args(3)%s, *) val
                    self%w(mult)%zero(ist1, ist1) = &
                    &    self%w(mult)%zero(ist1, ist1) + val
                end do
            case('kappa')
                call readf%next()
                call readf%parseline(' ')
                read(readf%args(1)%s, *) n_val
                do i = 1, n_val
                    call readf%next()
                    call readf%parseline(' ')
                    read(readf%args(1)%s, *) mult
                    read(readf%args(2)%s, *) ist1
                    read(readf%args(3)%s, *) imode1 
                    imode1 = imode1 - self%zero_mode
                    read(readf%args(4)%s, *) val
                    self%w(mult)%linear(imode1, ist1, ist1) = &
                    &    self%w(mult)%linear(imode1, ist1, ist1) + val
                end do
            case('lambda')
                call readf%next()
                call readf%parseline(' ')
                read(readf%args(1)%s, *) n_val
                do i = 1, n_val
                    call readf%next()
                    call readf%parseline(' ')
                    read(readf%args(1)%s, *) mult
                    read(readf%args(2)%s, *) ist1
                    read(readf%args(3)%s, *) ist2
                    read(readf%args(4)%s, *) imode1
                    imode1 = imode1 - self%zero_mode
                    read(readf%args(5)%s, *) val
                    self%w(mult)%linear(imode1, ist1, ist2) = &
                    &    self%w(mult)%linear(imode1, ist1, ist2) + val
                    if (ist1 == ist2) cycle
                    self%w(mult)%linear(imode1, ist2, ist1) = &
                    &    self%w(mult)%linear(imode1, ist2, ist1) + val
                end do
            case('gamma')
                call readf%next()
                call readf%parseline(' ')
                read(readf%args(1)%s, *) n_val
                do i = 1, n_val
                    call readf%next()
                    call readf%parseline(' ')
                    read(readf%args(1)%s, *) mult
                    read(readf%args(2)%s, *) ist1
                    read(readf%args(3)%s, *) imode1
                    imode1 = imode1 - self%zero_mode
                    read(readf%args(4)%s, *) imode2
                    imode2 = imode2 - self%zero_mode
                    read(readf%args(5)%s, *) val
                    self%w(mult)%quad(imode1, imode2, ist1, ist1) = &
                    &    self%w(mult)%quad(imode1, imode2, ist1, ist1) + val * 0.5_dp
                    if (imode1 == imode2) cycle
                    self%w(mult)%quad(imode2, imode1, ist1, ist1) = &
                    &    self%w(mult)%quad(imode2, imode1, ist1, ist1) + val * 0.5_dp
                end do
            case('SOC', 'DMX', 'DMY', 'DMZ')
                sec = readf%args(1)%s
                part = readf%args(2)%s
                if (part /= 'R') then
                    write(stderr, *) 'Error in vibronic_mod, vc_init subroutine.'
                    write(stderr, *) '  Only real SOC/DM matrices are supported.'
                    !> @todo Implement imaginary parts of SOC/DM matrices.
                    stop
                end if
                do i = 1, self%tot_ns
                    call readf%next()
                    read(readf%line, *) wrk(i, :)
                end do
                select case(sec)
                case('SOC')
                    call self%soc%add_zero(self%nmode, wrk)
                case('DMX')
                    call self%dm(1)%add_zero(self%nmode, wrk)
                case('DMY')
                    call self%dm(2)%add_zero(self%nmode, wrk)
                case('DMZ')
                    call self%dm(3)%add_zero(self%nmode, wrk)
                end select
            case default
                write(stderr, *) 'Error in vibronic_mod, vc_init subroutine.'
                write(stderr, *) '  Unrecognized keyword in template file: ', readf%args(1)%s
                stop
            end select            
        end do
        call readf%close()

        allocate(self%dw(self%nmode, 3))
        do i = 1, 3
            if (self%nstate(i) == 0) cycle
            do imode1 = 1, self%nmode
                call self%dw(imode1, i)%add_zero(self%nmode, self%w(i)%linear(imode1, :, :))
                call self%dw(imode1, i)%add_linear(2.0_dp * self%w(i)%quad(imode1, :, :, :))
            end do
        end do
    end subroutine vc_init


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_v0
    ! DESCRIPTION:
    !> @brief Read V0 file.
    !----------------------------------------------------------------------------------------------
    subroutine read_v0(self, v0_file)
        class(vc_model) :: self
        character(len=*), intent(in) :: v0_file
        type(reader) :: readf

        call readf%open(v0_file, abort_on_eof=.false.)

        do
            call readf%next()
            if (is_iostat_end(readf%iostat)) exit
            call readf%parseline(' ')
            select case(readf%args(1)%s)
            case('Frequencies')
                call readf%next()
                call readf%parseline(' ')
                self%nmode = readf%narg
                allocate(self%freq(self%nmode))
                read(readf%line, *) self%freq
            end select
        end do
        call readf%close()

        self%zero_mode = count(abs(self%freq) < 1.0e-12_dp)
        self%nmode = self%nmode - self%zero_mode
        if (self%zero_mode == 0) return
        if (self%zero_mode /= 6) then
            write(stderr, *) 'Warning in vibronic_mod, read_v0 subroutine.'
            write(stderr, *) '  Found ', self%zero_mode, ' zero frequencies. Expected none or 6.'
        end if
        
    end subroutine read_v0


end module vc_model_mod
