module timing_mod
    use global_defs
    implicit none

    private
    public :: timer

    !----------------------------------------------------------------------------------------------
    ! TYPE: timer
    !
    ! DESCRIPTION:
    !> @brief Object for measuring run time.
    !> @details 
    !! The start procedure sets the internal state of the timer. The print procedure outputs the
    !! time since the last call of start.
    !! Example use:
    !!
    !!     type(timer) :: clock
    !!     call clock%start()
    !!     ...
    !!     do_something
    !!     ...
    !!     call clock%print(stdout, message)
    !!
    !----------------------------------------------------------------------------------------------
    type timer
        integer(kind=8) :: t0 = 0 !< Start time.
    contains
        procedure :: start => timer_start !< Start the timer.
        procedure :: print => timer_print !< Print time since last start.
    end type timer


contains

    
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: timer_start
    !> @brief Initialize the timer.
    !----------------------------------------------------------------------------------------------
    subroutine timer_start(self)
        class(timer) :: self
        call system_clock(self%t0)
    end subroutine timer_start


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: timer_print
    !> @brief Print the current state of the timer.
    !> @details
    !! Time since the last update of the timer is printed in days/hours/minutes/seconds. The output
    !! is printed to outunit, preceeded by a message passed to the subroutine.
    !----------------------------------------------------------------------------------------------
    subroutine timer_print(self, outunit, msg)
        class(timer), intent(in) :: self
        integer, intent(in) :: outunit
        character(len=*), intent(in) :: msg
        integer(kind=8) :: now
        integer(kind=8) :: rate
        integer(kind=8) :: d, h, m
        real(dp) :: s
        call system_clock(now, rate)
        s = (now - self%t0) / real(rate, kind=dp)
        d = floor(s/86400.0_dp)
        s = s - d * 86400.0_dp
        h = floor(s/3600.0_dp)
        s = s - h * 3600.0_dp
        m = floor(s/60.0_dp)
        s = s - m * 60.0_dp
        write(outunit, 1000) msg, d, ' days, ', h, ' hours, ', m, ' minutes, ', s, ' seconds.'

1000 format (a,1x,i0,a,i0,a,i0,a,f5.2,a)

    end subroutine timer_print


end module timing_mod
