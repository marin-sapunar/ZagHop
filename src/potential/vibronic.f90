!--------------------------------------------------------------------------------------------------
! MODULE: vibronic_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date January, 2026
!
! DESCRIPTION: 
!> @brief Simple vibronic coupling interface.
!--------------------------------------------------------------------------------------------------
module vibronic_mod
    use global_defs
    use string_mod
    use file_mod, only : reader
    implicit none

    private
    public :: vc_model
    public :: vibronic_coupling

    
    type vc_property
        logical :: evaluated = .false.
        integer :: max_order = -1
        real(dp), allocatable :: zero(:, :)
        real(dp), allocatable :: linear(:, :, :)
        real(dp), allocatable :: quad(:, :, :, :)
        real(dp), allocatable :: prop(:, :)
        real(dp), allocatable :: trans_prop(:, :)
    contains
        procedure :: eval => vc_prop_eval
        procedure :: transform => vc_prop_transform
    end type vc_property


    type vc_model
        character(len=:), allocatable :: template_file
        character(len=:), allocatable :: v0_file
        integer :: nstate(3) = [0, 0, 0]
        integer :: nmode = 0
        integer :: tot_ns = 0
        integer, allocatable :: s2(:) !< 2x spin angular momentum quantuum numbers for each state.
        integer, allocatable :: ms2(:) !< 2x spin projection quantum numbers for each state.
        integer :: zero_mode = 0
        real(dp), allocatable :: freq(:)
        type(vc_property) :: h
        type(vc_property) :: soc
        type(vc_property) :: dm(3)
        type(vc_property), allocatable :: grad(:)
        real(dp), allocatable :: eigvec(:, :)
    contains
        procedure :: init => vc_init
        procedure :: read_v0
        procedure :: eval => evaluate_vc
        procedure :: get_grad => vc_get_grad
        procedure :: get_nadv => vc_get_nadv
    end type vc_model


    type(vc_model) :: vibronic_coupling


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
        integer :: n_val, i, j, k, cindex, mult
        integer ::  ist1, ist2, imode1, imode2
        integer :: i0(3) = [0, 0, 0]
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

        allocate(self%s2(self%tot_ns))
        allocate(self%ms2(self%tot_ns))
        cindex = 0
        do i = 1, 3
            i0(i) = cindex
            do j = 1, self%nstate(i)
                do k = -(i-1), (i-1), 2
                    cindex = cindex + 1
                    self%s2(cindex) = i - 1
                    self%ms2(cindex) = k
                end do
            end do
        end do

        self%h%max_order = 2
        allocate(self%h%zero(self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%h%linear(self%nmode, self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%h%quad(self%nmode, self%nmode, self%tot_ns, self%tot_ns), source=0.0_dp)

        allocate(wrk(self%tot_ns, self%tot_ns))
        do i = 1, self%nmode
            do j = 1, self%tot_ns
                self%h%quad(i, i, j, j) = 0.5_dp * self%freq(i + self%zero_mode)
            end do
        end do

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
                    ist1 = i0(mult) + mult * (ist1 - 1)
                    do j = 1, mult
                        self%h%zero(ist1 + j, ist1 + j) = val
                    end do
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
                    ist1 = i0(mult) + mult * (ist1 - 1)
                    do j = 1, mult
                        self%h%linear(imode1, ist1 + j, ist1 + j) = val
                    end do
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
                    ist1 = i0(mult) + mult * (ist1 - 1)
                    ist2 = i0(mult) + mult * (ist2 - 1)
                    do j = 1, mult
                        self%h%linear(imode1, ist1 + j, ist2 + j) = val
                        self%h%linear(imode1, ist2 + j, ist1 + j) = val
                    end do
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
                    ist1 = i0(mult) + mult * (ist1 - 1)
                    do j = 1, mult
                        self%h%quad(imode1, imode2, ist1+j, ist1+j) = &
                        &    self%h%quad(imode1, imode2, ist1+j, ist1+j) + val * 0.5_dp
                    end do
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
                    self%soc%max_order = 0
                    self%soc%zero = wrk
                    self%h%zero = self%h%zero + wrk
                case('DMX')
                    self%dm(1)%max_order = 0
                    self%dm(1)%zero = wrk
                case('DMY')
                    self%dm(2)%max_order = 0
                    self%dm(2)%zero = wrk
                case('DMZ')
                    self%dm(3)%max_order = 0
                    self%dm(3)%zero = wrk
                end select
            case default
                write(stderr, *) 'Error in vibronic_mod, vc_init subroutine.'
                write(stderr, *) '  Unrecognized keyword in template file: ', readf%args(1)%s
                stop
            end select            
        end do
        call readf%close()

        allocate(self%grad(self%nmode))
        do imode1 = 1, self%nmode
            self%grad(imode1)%max_order = self%h%max_order - 1
            allocate(self%grad(imode1)%zero, source=self%h%linear(imode1, :, :))
            allocate(self%grad(imode1)%linear, source=2.0_dp * self%h%quad(imode1, :, :, :))
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


    subroutine vc_prop_eval(self, q)
        class(vc_property) :: self
        real(dp), intent(in) :: q(:)
        integer :: i, j, imode, jmode

        if (self%evaluated) return

        if (.not. allocated(self%prop)) then
            allocate(self%prop(size(self%zero, 1), size(self%zero, 2)), source=self%zero)
        end if
        if (self%max_order == -1) then
            write(stderr, *) 'Error in vibronic_mod, vc_prop_eval subroutine.'
            write(stderr, *) '  Property not initialized.'
            stop
        end if

        self%prop = self%zero
        if (self%max_order >= 1) then
            do imode = 1, size(self%linear, 1)
                do i = 1, size(self%linear, 2)
                    do j = 1, size(self%linear, 3)
                        self%prop(i, j) = self%prop(i, j) + self%linear(imode, i, j) * q(imode)
                    end do
                end do
            end do
        end if
        if (self%max_order >= 2) then
            do imode = 1, size(self%quad, 1)
                do jmode = 1, size(self%quad, 2)
                    do i = 1, size(self%quad, 3)
                        do j = 1, size(self%quad, 4)
                            self%prop(i, j) = self%prop(i, j) + self%quad(imode, jmode, i, j) * q(imode) * q(jmode)
                        end do
                    end do
                end do
            end do
        end if
        self%evaluated = .true.
    end subroutine vc_prop_eval

    subroutine vc_prop_transform(self, trans_mat)
        class(vc_property) :: self
        real(dp), intent(in) :: trans_mat(:, :)

        if (.not. allocated(self%prop)) then
            write(stderr, *) 'Error in vibronic_mod, vc_prop_transform subroutine.'
            write(stderr, *) '  Property not evaluated yet.'
            stop
        end if

        self%trans_prop = matmul(transpose(trans_mat), matmul(self%prop, trans_mat))
    end subroutine vc_prop_transform


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: evaluate_vc
    ! DESCRIPTION:
    !> @brief Evaluate vibronic coupling Hamiltonian at given nuclear geometry.
    !----------------------------------------------------------------------------------------------
    subroutine evaluate_vc(self, q, basis, adiab_e)
        use linalg_wrapper_mod, only : syev
        class(vc_model) :: self
        real(dp), intent(in) :: q(:)
        character(len=*), intent(in) :: basis
        real(dp), allocatable, intent(out) :: adiab_e(:)
        integer :: i, j, i0, imode
        real(dp), allocatable :: wrk(:, :), wrk_e(:)

        self%h%evaluated = .false.
        do imode = 1, self%nmode
            self%grad(imode)%evaluated = .false.
        end do
        if (self%soc%max_order >= 0) then
            self%soc%evaluated = .false.
        end if
        do i = 1, 3
            if (self%dm(i)%max_order >= 0) then
                self%dm(i)%evaluated = .false.
            end if
        end do

        if (.not. allocated(adiab_e)) then
            allocate(adiab_e(self%tot_ns), source=0.0_dp)
        end if
        if (.not. allocated(self%eigvec)) then
            allocate(self%eigvec(self%tot_ns, self%tot_ns), source=0.0_dp)
        end if

        call self%h%eval(q)
        do imode = 1, self%nmode
            call self%grad(imode)%eval(q)
        end do
        if (self%soc%max_order >= 0) then
            call self%soc%eval(q)
        end if
        do i = 1, 3
            if (self%dm(i)%max_order >= 0) then
                call self%dm(i)%eval(q)
            end if
        end do

        wrk = self%h%prop
        select case(basis)
        case('adiabatic')
            call syev(wrk, adiab_e, jobz='V', uplo='U')
            self%eigvec = wrk
        case('spin-diabatic')
            i0 = 0
            do i = 1, 3
                if (self%nstate(i) == 0) cycle
                if (allocated(wrk)) deallocate(wrk)
                allocate(wrk(self%nstate(i), self%nstate(i)))
                if (allocated(wrk_e)) deallocate(wrk_e)
                allocate(wrk_e(self%nstate(i)))
                wrk = self%h%prop(i0+1:i0+i*self%nstate(i):i, i0+1:i0+i*self%nstate(i):i)
                call syev(wrk, wrk_e, jobz='V', uplo='U')
                do j = 1, i
                    adiab_e(i0+j:i0+i*self%nstate(i):i) = wrk_e
                    self%eigvec(i0+j:i0+i*self%nstate(i):i, i0+j:i0+i*self%nstate(i):i) = wrk
                end do
                i0 = i0 + i * self%nstate(i)
            end do
        end select

        call self%h%transform(self%eigvec)
        do imode = 1, self%nmode
            call self%grad(imode)%transform(self%eigvec)
        end do
        if (self%soc%max_order >= 0) then
            call self%soc%transform(self%eigvec)
        end if
        do i = 1, 3
            if (self%dm(i)%max_order >= 0) then
                call self%dm(i)%transform(self%eigvec)
            end if
        end do

    end subroutine evaluate_vc

    function vc_get_grad(self, istate) result(grad)
        class(vc_model) :: self
        integer, intent(in) :: istate
        real(dp), allocatable :: grad(:)
        integer :: imode

        allocate(grad(self%nmode), source=0.0_dp)
        do imode = 1, self%nmode
            grad(imode) = self%grad(imode)%trans_prop(istate, istate)
        end do
    end function vc_get_grad

    function vc_get_nadv(self, istate1, istate2) result(nadv)
        class(vc_model) :: self
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)
        integer :: imode
        real(dp) :: edif
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        allocate(nadv(self%nmode), source=0.0_dp)
        edif = self%h%trans_prop(istate2, istate2) - self%h%trans_prop(istate1, istate1)
        if (abs(edif) < tiny_hf) then
            if (stdp1) then
                write(stderr, *) 'Warning in vibronic_mod, evaluate_vc subroutine.'
                write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
            end if
            edif = sign(tiny_hf, edif)
        end if
        do imode = 1, self%nmode
            nadv(imode) = self%grad(imode)%trans_prop(istate1, istate2) / edif
        end do

    end function vc_get_nadv


end module vibronic_mod
