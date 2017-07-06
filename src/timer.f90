!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module timer
  implicit none
  private

  integer :: itime_ini
  logical :: if_initialized = .false.

  public :: init_timer, &
            elapse_time

contains
!-----------------------------------------------------------------------------------------
  subroutine init_timer
    implicit none

    call system_clock(itime_ini)
    if_initialized = .true.

  end subroutine init_timer
!-----------------------------------------------------------------------------------------
  subroutine elapse_time(time)
    implicit none
    real(8),intent(out) :: time
    integer :: itime,itime_rate,itime_max

    if(.not.if_initialized)write(*,"(A)") &
      "Warning: time is not initialized. The elapsed time can be wrong."

    call system_clock(itime,itime_rate,itime_max)

    if(itime < itime_ini)then
      time = dble(itime_max-itime_ini + itime + 1)/itime_rate
    else
      time = dble(itime - itime_ini)/itime_rate
    end if
    

  end subroutine elapse_time
!-----------------------------------------------------------------------------------------
end module timer
