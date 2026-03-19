program run_vc
    use global_defs
    use vibronic_mod, only : vc_model
    implicit none

    type(vc_model) :: vc
    character(len=256) :: template_file, geom_file
    real(dp), allocatable :: q(:), energies(:)
    integer :: i, nmode, io, iunit
    character(len=32) :: label, tag
    real(dp) :: mass, qi
    real(dp), allocatable :: oscill(:)

    if (command_argument_count() < 2) then
        write(stderr, '(a)') 'Usage: run_vc <template_file> <geom_file>'
        stop 1
    end if

    call get_command_argument(1, template_file)
    call get_command_argument(2, geom_file)

    call vc%init(trim(template_file))
    nmode = vc%nmode

    allocate(q(nmode))

    open(newunit=iunit, file=trim(geom_file), status='old', action='read', iostat=io)
    if (io /= 0) then
        write(stderr, '(a)') 'Error: cannot open geom file: '//trim(geom_file)
        stop 1
    end if
    do i = 1, nmode
        read(iunit, *, iostat=io) label, mass, q(i), tag
        if (io /= 0) then
            write(stderr, '(a,i0,a)') 'Error reading mode ', i, ' from geom file.'
            stop 1
        end if
    end do
    close(iunit)

    call vc%eval(q, 'spin-diabatic', energies)

    if (vc%dm(1)%max_order >= 0) then
        oscill = vc%get_oscill(1)
    end if
    open(newunit=iunit, file='qm_en.dat', status='replace', action='write', iostat=io)
    do i = 1, vc%tot_ns
        if (vc%dm(1)%max_order >= 0) then
            write(iunit, '(i4, f20.12, f20.12)') i, energies(i), oscill(i)
        else
            write(iunit, '(i4, f20.12)') i, energies(i)
        end if
    end do
    close(iunit)

    open(newunit=iunit, file='qm_adt.dat', status='replace', action='write', iostat=io)
    do i = 1, vc%tot_ns
        write(iunit, '(*(e20.12))') vc%eigvec(:, i)
    end do
    close(iunit)

end program run_vc
