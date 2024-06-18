!--------------------------------------------------------------------------------------------------
! MODULE: StringMod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Helper routines for working with strings.
!--------------------------------------------------------------------------------------------------
module string_mod
    use global_defs
    implicit none
    
    private
    public :: string
    public :: string_parse
    public :: toupper
    public :: tolower
    public :: read_index_list    
    public :: read_index_list_unsort
    public :: char_is_num
    public :: rmwhite
    public :: compact

    
    type string
        character(len=:), allocatable :: s
    end type string


    character(len=10), parameter :: nums = '1234567890'

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: string_parse
    !
    ! DESCRIPTION:
    !> @brief 
    !! Parse a string into an array of arguments args(1), ..., args(nargs).
    !> @details 
    !! The string is separated based on the delimiters contained in string 'delims'. Consecutive
    !! delimiters are treated as one. Whitespace at the beginning and end of each argument is 
    !! removed.
    !----------------------------------------------------------------------------------------------
    subroutine string_parse(input_string, delimiters, narg, arguments, protect_quotes)
        character(len=*), intent(in) :: input_string
        character(len=*), intent(in) :: delimiters !< String containing the possible delimiters between arguments.
        integer, intent(out) :: narg
        type(string), allocatable, intent(out) :: arguments(:)
        logical, intent(in), optional :: protect_quotes
      
        character(len=:), allocatable :: tempstr
        character(len=:), allocatable :: tempout
        integer :: i
        logical :: pquotes
      
      
        tempstr = input_string
        call compact(tempstr)
        pquotes = .true.
        IF (present(protect_quotes)) pquotes = protect_quotes
        narg = 0
      
        ! First determine number of arguments.
        do
            if (len_trim(tempstr) == 0) exit
            narg = narg + 1
            call split(tempstr, delimiters, tempout, pquotes)
        end do
        ! Allocate and fill arguments array.
        if (allocated(arguments)) deallocate(arguments)
        allocate(arguments(narg))
        tempstr = input_string
        call compact(tempstr)
        do i = 1, narg
            call split(tempstr, delimiters, tempout, pquotes)
            arguments(i)%s = tempout
        end do
    end subroutine string_parse
 
 
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_index_list
    !
    ! DESCRIPTION:
    !> @brief Read a list of positive integers from a string.
    !> @details 
    !! The integers can be given as a comma/whitespace separated list. Members of the list can be 
    !! single integers or ranges separated by a hyphen.
    !----------------------------------------------------------------------------------------------
    subroutine read_index_list(str, indexlist)
        use sort_mod, only : sort
        character(len=*), intent(in) :: str !< Input string.
        integer, allocatable, intent(out) :: indexlist(:) !< Final integer list.
       
        character(len=:), allocatable :: tempstr
        logical :: range_0, range_1
        integer :: i, j, i0, last
        integer :: ntot
        integer :: int1, int0
        integer :: templist(20000)
       
        ntot = 0
        tempstr = str
        call compact(tempstr)

        if (len(tempstr) == 1) then
            allocate(indexlist(1))
            read(tempstr, *) indexlist(1)
            return
        end if
       
        i0 = 1
        range_0 = .false.
        range_1 = .false.
        do i = 2, len(tempstr)
            select case(tempstr(i:i))
            case('-')
                last = 1
                range_0 = .true.
            case(',')
                last = 2
            case(' ')
                if (last /= 4) cycle
            case default
                last = 4
                if (i /= len(tempstr)) cycle
            end select
            if (i /= len(tempstr)) then ! Reached a delimiter.
                read(tempstr(i0:i-1), *) int1
            else ! Reached end of string.
                read(tempstr(i0:i), *) int1
            end if
            i0 = i + 1
            if (range_0) then
                int0 = int1
                range_0 = .false.
                range_1 = .true.
            else if (range_1) then
                range_1 = .false.
                do j = int0, int1
                    ntot = ntot + 1
                    templist(ntot) = j
                end do
            else
                ntot = ntot + 1
                templist(ntot) = int1
            end if
        end do

       
        if (allocated(indexlist)) deallocate(indexlist)
        allocate(indexlist(ntot))
        call sort(templist(1:ntot))
        indexlist = templist(1:ntot)
 
    end subroutine read_index_list
 
     !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_index_list
    !
    ! DESCRIPTION:
    !> @brief Read a list of positive integers from a string.
    !> @details 
    !! The integers can be given as a comma/whitespace separated list. Members of the list can be 
    !! single integers or ranges separated by a hyphen.
    !----------------------------------------------------------------------------------------------
    subroutine read_index_list_unsort(str, indexlist)
        use sort_mod, only : sort
        character(len=*), intent(in) :: str !< Input string.
        integer, allocatable, intent(out) :: indexlist(:) !< Final integer list.
       
        character(len=:), allocatable :: tempstr
        logical :: range_0, range_1
        integer :: i, j, i0, last
        integer :: ntot
        integer :: int1, int0
        integer :: templist(20000)
       
        ntot = 0
        tempstr = str
        call compact(tempstr)

        if (len(tempstr) == 1) then
            allocate(indexlist(1))
            read(tempstr, *) indexlist(1)
            return
        end if
       
        i0 = 1
        range_0 = .false.
        range_1 = .false.
        do i = 2, len(tempstr)
            select case(tempstr(i:i))
            case('-')
                last = 1
                range_0 = .true.
            case(',')
                last = 2
            case(' ')
                if (last /= 4) cycle
            case default
                last = 4
                if (i /= len(tempstr)) cycle
            end select
            if (i /= len(tempstr)) then ! Reached a delimiter.
                read(tempstr(i0:i-1), *) int1
            else ! Reached end of string.
                read(tempstr(i0:i), *) int1
            end if
            i0 = i + 1
            if (range_0) then
                int0 = int1
                range_0 = .false.
                range_1 = .true.
            else if (range_1) then
                range_1 = .false.
                do j = int0, int1
                    ntot = ntot + 1
                    templist(ntot) = j
                end do
            else
                ntot = ntot + 1
                templist(ntot) = int1
            end if
        end do

       
        if (allocated(indexlist)) deallocate(indexlist)
        allocate(indexlist(ntot))
!        call sort(templist(1:ntot))
        indexlist = templist(1:ntot)
 
    end subroutine read_index_list_unsort
 
 
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rmWhiteSpace
    !
    ! DESCRIPTION:
    !> @brief 
    !! Remove all spaces, tabs and control characters from a string.
    !----------------------------------------------------------------------------------------------
    SUBROUTINE rmWhite(str)
       character(len=:), intent(inout), allocatable :: str
 
       character(len=:), allocatable :: outstr
       character(len=1) :: ch
       integer :: i
       integer :: ich
       str = adjustl(str)
       outstr = ''
 
       DO i = 1, len_trim(str)
          ch = str(i:i)
          ich = iachar(ch)
          SELECT CASE (ich)
             case(33:)
                outstr = outstr // ch
          END SELECT
       END DO
 
       str = adjustl(outstr)
 
    END SUBROUTINE rmWhite
 

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: compact
    !
    ! DESCRIPTION:
    !> @brief 
    !! Convert multiple spaces and tabs to single spaces and remove control characters from a 
    !! string.
    !----------------------------------------------------------------------------------------------
    SUBROUTINE compact(str)
        character(len=:), intent(inout), allocatable :: str
      
        character(len=:), allocatable :: outstr
        character(len=1) :: ch
        integer :: i
        integer :: isp
        integer :: ich
      
        str = adjustl(str)
        outstr = ''
        isp = 1
        
        DO i = 1, len_trim(str)
            ch = str(i:i)
            ich = iachar(ch)
            SELECT CASE (ich)
                case(9,11,32)
                    IF (isp == 0) outstr = outstr // ' '
                    isp = 1
                case(33:)
                    outstr = outstr // ch
                    isp = 0
            END SELECT
        END DO
      
        str = adjustl(outstr)
        
    END SUBROUTINE compact
 
    
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: split
    !
    ! DESCRIPTION:
    !> @brief 
    !! Search input string for first occurence of any character from delimiters string and split
    !! it at the position.
    !----------------------------------------------------------------------------------------------
    SUBROUTINE split(str, delims, before, protect_quotes)
        character(len=:), intent(inout), allocatable :: str !< Input string. Second part of output string.
        character(len=*), intent(in) :: delims !< Characters at which the string is split.
        character(len=:), intent(out), allocatable :: before !< First part of output string.
        logical, intent(in) :: protect_quotes
        
        character :: ch
        integer :: i
        integer :: ipos
        logical :: found, quoted
        character(len=:), allocatable :: tmpstr
        
        found = .false.
        str = adjustl(str)
        before = ''
        quoted = .false.
        DO i = 1, len_trim(str)
            ch = str(i:i)
            IF ((i == 1) .AND. protect_quotes) THEN
                IF ((ch == '"') .OR. (ch == "'")) THEN
                    quoted = .true.
                    CYCLE
                END IF
            END IF
            IF (quoted) THEN
                IF ((ch == '"') .OR. (ch == "'")) THEN
                    quoted = .false.
                    CYCLE
                END IF
            END IF
            ipos = index(delims, ch)
            IF ((ipos == 0) .OR. quoted) THEN
                IF (found) THEN
                    ipos = len(str)
                    tmpstr = str(i:len(str))
                    str = tmpstr
                    EXIT
                END IF
                before = before // ch
            ELSE
                IF (before == '') CYCLE
                IF (i == len_trim(str)) THEN
                    found = .false.
                    EXIT
                END IF
                found = .true.
            END IF
        END DO
        
        IF (.not. found) str = ''
        
    END SUBROUTINE split

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ToUpper
    !
    ! DESCRIPTION:
    !> @brief 
    !! Convert all letters in a string to upper case. (Only works with ASCII characters).
    !----------------------------------------------------------------------------------------------
    function ToUpper(strIn) result(strOut)
        character(len=*), intent(in) :: strIn
        character(len=:), allocatable :: strOut
        integer :: i
        integer :: j
        integer :: offset
  
        offset=iachar('A')-iachar('a')
        strOut = strIn
        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j>= iachar("a") .and. j<=iachar("z") ) then
                strOut(i:i) = achar(iachar(strIn(i:i)) + offset)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do
    end function ToUpper

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ToLower
    !
    ! DESCRIPTION:
    !> @brief 
    !! Convert all letters in a string to lower case. (Only works with ASCII characters).
    !----------------------------------------------------------------------------------------------
    pure function ToLower(strIn) result(strOut)
        character(len=*), intent(in) :: strIn
        character(len=:), allocatable :: strOut
        integer :: i
        integer :: j
        integer :: offset
 
        offset=iachar('A')-iachar('a')
        strOut = strIn
        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j>= iachar("A") .and. j<=iachar("Z") ) then
                strOut(i:i) = achar(iachar(strIn(i:i)) - offset)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do
    end function ToLower


    pure function char_is_num(cha) result(is_numeric)
        character(len=*), intent(in) :: cha
        logical :: is_numeric

        is_numeric = (verify(cha, nums) == 0)

    end function char_is_num


end module string_mod
