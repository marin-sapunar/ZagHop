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
        real(dp), allocatable :: zero_order(:, :)
        real(dp), allocatable :: linear(:, :, :)
        real(dp), allocatable :: quadratic(:, :, :, :)
        real(dp), allocatable :: soc(:, :)
        real(dp), allocatable :: dm(:, :, :)
        real(dp), allocatable :: diab_h(:, :)
        real(dp), allocatable :: eigvec(:, :)
        real(dp), allocatable :: diab_grad(:, :, :)
        real(dp), allocatable :: adiab_grad(:, :, :)
    contains
        procedure :: init => vc_init
        procedure :: read_v0
        procedure :: eval => evaluate_vc
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

        allocate(wrk(self%tot_ns, self%tot_ns))
        allocate(self%zero_order(self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%linear(self%nmode, self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%soc(self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%dm(3, self%tot_ns, self%tot_ns), source=0.0_dp)
        allocate(self%quadratic(self%nmode, self%nmode, self%tot_ns, self%tot_ns), source=0.0_dp)
        do i = 1, self%nmode
            do j = 1, self%tot_ns
                self%quadratic(i, i, j, j) = 0.5_dp * self%freq(i + self%zero_mode)
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
                        self%zero_order(ist1 + j, ist1 + j) = val
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
                        self%linear(imode1, ist1 + j, ist1 + j) = val
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
                        self%linear(imode1, ist1 + j, ist2 + j) = val
                        self%linear(imode1, ist2 + j, ist1 + j) = val
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
                        self%quadratic(imode1, imode2, ist1+j, ist1+j) = &
                        &    self%quadratic(imode1, imode2, ist1+j, ist1+j) + val * 0.5_dp
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
                    self%soc = wrk
                case('DMX')
                    self%dm(1, :, :) = wrk
                case('DMY')
                    self%dm(2, :, :) = wrk
                case('DMZ')
                    self%dm(3, :, :) = wrk
                end select
            case default
                write(stderr, *) 'Error in vibronic_mod, vc_init subroutine.'
                write(stderr, *) '  Unrecognized keyword in template file: ', readf%args(1)%s
                stop
            end select            
        end do
        call readf%close()
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
        integer :: i, j, imode, jmode
        integer :: i0, i_end
        real(dp), allocatable :: wrk(:, :), wrk_e(:)
        real(dp) :: edif
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        if (.not. allocated(adiab_e)) then
            allocate(adiab_e(self%tot_ns), source=0.0_dp)
        end if
        if (.not. allocated(self%adiab_grad)) then
            allocate(self%adiab_grad(self%nmode, self%tot_ns, self%tot_ns), source=0.0_dp)
        end if
        if (.not. allocated(self%eigvec)) then
            allocate(self%eigvec(self%tot_ns, self%tot_ns), source=0.0_dp)
        end if

        wrk = self%zero_order
        wrk = wrk + self%soc

        do imode = 1, self%nmode
            do i = 1, self%tot_ns
                do j = 1, self%tot_ns
                    wrk(i, j) = wrk(i, j) + self%linear(imode, i, j) * q(imode)
                end do
            end do
            do jmode = 1, self%nmode
                do i = 1, self%tot_ns
                    do j = 1, self%tot_ns
                        wrk(i, j) = wrk(i, j) + self%quadratic(imode, jmode, i, j) * q(imode) * q(jmode)
                    end do
                end do
            end do
        end do


        self%diab_h = wrk
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
                wrk = self%diab_h(i0+1:i0+i*self%nstate(i):i, i0+1:i0+i*self%nstate(i):i)
                call syev(wrk, wrk_e, jobz='V', uplo='U')
                do j = 1, i
                    adiab_e(i0+j:i0+i*self%nstate(i):i) = wrk_e
                    self%eigvec(i0+j:i0+i*self%nstate(i):i, i0+j:i0+i*self%nstate(i):i) = wrk
                end do
                i0 = i0 + i * self%nstate(i)
            end do
        end select


        self%diab_grad = self%linear
        do imode = 1, self%nmode
            do i = 1, self%tot_ns
                do j = 1, self%tot_ns
                    do jmode = 1, self%nmode
                        self%diab_grad(imode, i, j) = self%diab_grad(imode, i, j) + &
                            2.0_dp * self%quadratic(imode, jmode, i, j) * q(jmode)
                    end do
                end do
            end do
        end do

        do imode = 1, self%nmode
            self%adiab_grad(imode, :, :) = matmul(matmul(transpose(self%eigvec), &
            &                                     self%diab_grad(imode, :, :)), self%eigvec)
            do i = 1, self%tot_ns
                do j = 1, self%tot_ns
                    if (i==j) cycle
                    if ((self%s2(i) /= self%s2(j)) .or. (self%ms2(i) /= self%ms2(j))) cycle
                    edif = adiab_e(j) - adiab_e(i)
                    if (abs(edif) < tiny_hf) then
                        if (stdp1) then
                            write(stderr, *) 'Warning in vibronic_mod, evaluate_vc subroutine.'
                            write(stderr, *) '  Near-degeneracy between states ', i, ' and ', j, '.'
                            write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
                        end if
                        edif = sign(tiny_hf, edif)
                    end if
                    self%adiab_grad(imode, i, j) = self%adiab_grad(imode, i, j) / edif
                end do
            end do
        end do
    end subroutine evaluate_vc

end module vibronic_mod
