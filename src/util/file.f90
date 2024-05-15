!--------------------------------------------------------------------------------------------------
! MODULE: file_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Helper subroutines for reading and parsing lines of text from input files.
!--------------------------------------------------------------------------------------------------
module file_mod
    use global_defs
    use string_mod
    implicit none
    
    private
    public :: reader
    public :: check_string_in_file
    public :: check_is_dir
    public :: need_dir
    public :: need_file
    

    !----------------------------------------------------------------------------------------------
    ! TYPE: Reader
    !
    !> @brief Class for reading and parsing files.
    !> @details 
    !! Class to open a file and read from it. Contains procedure (%next) to read arbitrary length
    !! lines while ignoring comments (and empty lines), and procedure (%parse) to split current
    !! line into arguments based on delimiters. Also contains subroutines for searching within the
    !! file.
    !----------------------------------------------------------------------------------------------
    type reader
        character(len=:), allocatable :: file
        integer :: unit
        character(len=:), allocatable :: comment
        character(len=:), allocatable :: continuation
        character(len=:), allocatable :: delimiters
        logical :: skip_empty = .true.
        logical :: abort_on_eof = .true.
        logical :: abort_on_error = .true.
        integer :: iostat = 0
        integer :: line_num = 0
        character(len=:), allocatable :: line
        type(string), allocatable :: args(:)
        integer :: narg
    contains
        procedure :: open => reader_open
        procedure :: close => reader_close
        procedure :: rewind => reader_rewind
        procedure :: next => reader_next
        procedure :: parseline => reader_parseline
        procedure :: go_to_line => reader_go_to_line
        procedure :: go_to_keyword => reader_go_to_keyword
        procedure :: count_lines_to_keyword => reader_count_lines_to_keyword
        procedure :: count_keyword_appearances => reader_count_keyword_appearances
    end type reader


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_open
    !
    !> @brief Initialize an instance of the reader class.
    !> @details
    !! Open a file within a reader class. Also sets options used by other subroutines when reading 
    !! the file.
    !----------------------------------------------------------------------------------------------
    subroutine reader_open(self, file_name, comment, delimiters, continuation, skip_empty,        &
                           abort_on_eof, abort_on_error)
        class(reader) :: self !< Reader instance to initialize.
        character(len=*), intent(in) :: file_name !< Name of file.
        character(len=*), intent(in), optional :: comment !< Characters starting a comment to skip.
        character(len=*), intent(in), optional :: delimiters !< Delimiter characters.
        character(len=*), intent(in), optional :: continuation !< Line continuation characters.
        logical, intent(in), optional :: skip_empty !< Skip empty lines when reading.
        logical, intent(in), optional :: abort_on_eof !< Abort on eof or return iostat.
        logical, intent(in), optional :: abort_on_error !< Abort on eof or return iostat.

        self%file=trim(file_name)
        call need_file(self%file)
        open(newunit=self%unit, file=self%file, action='read', status='old')
    
        if (.not. allocated(self%comment)) self%comment = ''
        if (.not. allocated(self%continuation)) self%continuation = ''
        if (.not. allocated(self%delimiters)) self%delimiters = ' '
        if (present(comment)) self%comment = comment
        if (present(delimiters)) self%delimiters = delimiters
        if (present(continuation)) self%continuation = continuation
        if (present(skip_empty)) self%skip_empty = skip_empty
        if (present(abort_on_eof)) self%abort_on_eof = abort_on_eof
        if (present(abort_on_error)) self%abort_on_error = abort_on_error
    end subroutine reader_open


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_close
    !> @brief Close opened reader instance.
    !----------------------------------------------------------------------------------------------
    subroutine reader_close(self)
        class(reader) :: self

        close(self%unit)
        self%file = ''
        self%iostat = 0
        self%line_num = 0
    end subroutine reader_close


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_rewind
    !> @brief Rewind opened reader instance.
    !----------------------------------------------------------------------------------------------
    subroutine reader_rewind(self)
        class(reader) :: self

        rewind(self%unit)
        self%iostat = 0
        self%line_num = 0
    end subroutine reader_rewind


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_next
    !
    !> @brief Read current line from an open unit.
    !> @details 
    !! Can handle lines of arbitrary length. Uses options set in the reader class instance.
    !----------------------------------------------------------------------------------------------
    subroutine reader_next(self, abort_on_error, abort_on_eof)
        class(reader) :: self !< Reader instance.
        logical, intent(in), optional :: abort_on_error !< Abort on eof or return iostat.
        logical, intent(in), optional :: abort_on_eof !< Abort on eof or return iostat.
      
        integer :: i
        integer :: sz
        integer :: ipos
        character(len=256) :: buffer
        logical :: abort_error
        logical :: abort_eof

      
        if (present(abort_on_error)) then
            abort_error = abort_on_error
        else
            abort_error = self%abort_on_error
        end if
        if (present(abort_on_eof)) then
            abort_eof = abort_on_eof
        else
            abort_eof = self%abort_on_eof
        end if

        self%line = ''
101     self%line_num = self%line_num + 1
        do
            read(self%unit, "(a)", advance='no', iostat=self%iostat, size=sz) buffer
            if (self%iostat == 0) then
                self%line = self%line // buffer(:sz)
            else if (is_iostat_eor(self%iostat)) then
                self%line = self%line // buffer(:sz)
                exit
            else if (is_iostat_end(self%iostat)) then
                if (abort_eof) then
                    write(stderr, *) 'Error in read_line subroutine. Unexpected end of file.'
                    call abort()
                end if
                return
            else
                if (abort_error) then
                    write(stderr, '(1x,a,i0)') 'Error in read_line subroutine. Iostat: ', self%iostat
                    call abort()
                end if
                return
            end if
        end do
      
        ! Remove comments from line.
        do i = 1, len_trim(self%comment)
            ipos = index(self%line, self%comment(i:i))
            if (ipos /= 0) self%line = self%line(:ipos-1)
        end do
        ! Check if line is empty.
        if (self%skip_empty) then
            if (self%line == '') goto 101
        end if
        ! If line ends with a continuation character, continue reading next line.
        ipos = len_trim(self%line)
        do i = 1, len_trim(self%continuation)
            if (self%line(ipos:ipos) == self%continuation(i:i)) then
                goto 101
            end if
        end do
    end subroutine reader_next


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_parseline
    !
    !> @brief Parse the current line in a reader variable.
    !> @details 
    !! The line is separated into strings based on the delimiters. Consecutive delimiters are
    !! treated as one. Whitespace at the beginning and end of each argument is removed.
    !----------------------------------------------------------------------------------------------
    subroutine reader_parseline(self, delimiters)
        use string_mod, only : string_parse
        class(reader) :: self !< Reader instance.
        character(len=*), intent(in), optional :: delimiters !< Delimiters separating arguments.
      
        if (present(delimiters)) then
            call string_parse(self%line, delimiters, self%narg, self%args)
        else
            call string_parse(self%line, self%delimiters, self%narg, self%args)
        end if
    end subroutine reader_parseline


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_go_to_keyword
    !
    !> @brief Go to the next line containing a keyword.
    !> @details
    !! Go through a file using the reader class and stop when a line containing the given keyword is
    !! found. The reader can then be used to read through this section. If the keyword is an empty
    !! string, the subroutine finds the next empty line.
    !----------------------------------------------------------------------------------------------
    subroutine reader_go_to_keyword(self, keyword, found, count_lines, case_sensitive)
        class(reader) :: self !< Reader instance.
        character(len=*), intent(in) :: keyword !< Keyword to search for.
        logical, intent(out), optional :: found !< Status on exit. .true. if keyword is found, and 
            !! .false. otherwise. If not present, the program is aborted if keyword is not found.
        integer, intent(out), optional :: count_lines !< Number of lines skipped during search.
        logical, intent(in), optional :: case_sensitive !< Should the search be case sensitive.
        logical :: tcase_sensitive
        logical :: check
        integer :: c
        character(len=:), allocatable :: errmsg
        
        tcase_sensitive = .true.
        if (present(case_sensitive)) tcase_sensitive = case_sensitive

        c = -1
        do
            call self%next(abort_on_eof = .false.)
            c = c + 1
            if (is_iostat_end(self%iostat)) then
                if (present(found)) then
                    found = .false.
                    return
                else
                    errmsg = 'Keyword "'//keyword//'" not found in file "'//self%file//'".'
                    call errstop('reader_go_to_keyword', errmsg, 1)
                end if
            end if
            if (keyword == '') then
                check = (self%line == keyword)
            else if (tcase_sensitive) then
                check = (index(self%line, keyword) /= 0)
            else
                check = (index(tolower(self%line), tolower(keyword)) /= 0)
            end if
            if (check) then
                if (present(found)) found = .true.
                if (present(count_lines)) count_lines = c
                return
            end if
        end do
    end subroutine reader_go_to_keyword


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: reader_go_to_line
    !> @brief Go to line number in reader file.
    !----------------------------------------------------------------------------------------------
    subroutine reader_go_to_line(self, line)
        class(reader) :: self
        integer, intent(in) :: line
        integer :: i

        if (self%line_num == line) then
            return
        else if (self%line_num > line) then
            call self%rewind()
        end if
        do i = self%line_num, line
            call self%next()
            if (self%line_num == line) exit
        end do
    end subroutine reader_go_to_line


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: reader_count_lines_to_keyword
    !> @brief Count number of lines remaining until next appearence of given keyword.
    !----------------------------------------------------------------------------------------------
    function reader_count_lines_to_keyword(self, keyword) result(nlines)
        class(reader) :: self !< Reader instance.
        character(len=*), intent(in) :: keyword !< Keyword to search for.
        integer :: nlines
        integer :: iline

        iline = self%line_num
        call self%go_to_keyword(keyword, count_lines=nlines)
        call self%go_to_line(iline)
    end function reader_count_lines_to_keyword


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: reader_count_keyword_appearances
    !> @brief Count number of appearances of keyword until end of file.
    !----------------------------------------------------------------------------------------------
    function reader_count_keyword_appearances(self, keyword) result(nkey)
        class(reader) :: self !< Reader instance.
        character(len=*), intent(in) :: keyword !< Keyword to search for.
        integer :: nkey
        integer :: iline
        logical :: chk

        iline = self%line_num
        nkey = 0
        do
            call self%go_to_keyword(keyword, found=chk)
            if (.not. chk) exit
            nkey = nkey + 1
        end do
        call self%go_to_line(iline)
    end function reader_count_keyword_appearances


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_string_in_file
    !> @brief Return .true. if string is present anywhere in file. Return .false. otherwise.
    !----------------------------------------------------------------------------------------------
    function check_string_in_file(file_name, keyword) result(flag)
        character(len=*), intent(in) :: file_name
        character(len=*), intent(in) :: keyword
        logical :: flag
        type(reader) :: reader_file


        flag = .false.
        call need_file(file_name, 'While checking for string: '//keyword)
        call reader_file%open(file_name)
        call reader_file%go_to_keyword(keyword, found=flag)
        call reader_file%close
    end function check_string_in_file


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_is_dir
    !> @brief Return .true. if given string is a valid path to a directory.
    !----------------------------------------------------------------------------------------------
    function check_is_dir(test_dir) result(is_dir)
        character(len=*), intent(in) :: test_dir
        logical :: is_dir
#if __INTEL_COMPILER
        inquire(directory=test_dir, exist=is_dir)
#else
        inquire(file=trim(test_dir)//'/.', exist=is_dir)
#endif
    end function check_is_dir


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: need_file
    !> @brief Print error message and abort program if file is not found.
    !----------------------------------------------------------------------------------------------
    subroutine need_file(test_file, error_message)
        character(len=*), intent(in) :: test_file
        character(len=*), intent(in), optional :: error_message
        logical :: chk

        inquire(file=test_file, exist=chk)
        if (.not. chk) then
            write(stderr, '(1x,a,a,a)') 'Error. File ', test_file, ' not found.'
            if (present(error_message)) write(stderr, *) error_message
            stop
        end if
    end subroutine need_file


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: need_dir
    !> @brief Print error message and abort program if directory is not found.
    !----------------------------------------------------------------------------------------------
    subroutine need_dir(test_dir, error_message)
        character(len=*), intent(in) :: test_dir
        character(len=*), intent(in), optional :: error_message

        if (.not. check_is_dir(test_dir)) then
            write(stderr, *)
            write(stderr, '(1x,a,a,a)') 'Error. Directory ', test_dir, ' not found.'
            if (present(error_message)) write(stderr, *) error_message
            stop
        end if
    end subroutine need_dir


end module file_mod
