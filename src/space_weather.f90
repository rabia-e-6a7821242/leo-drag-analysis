! space_weather.f90
! Solar (F10.7) and geomagnetic (Ap) drivers as a function of time.
!
! solar_mode = 'constant' -> F10.7 = f107_const for the whole run
! solar_mode = 'cycle'    -> idealised 11-year cycle:
!     F10.7(t) = f107_min + (f107_max - f107_min) * sin^2( pi (t + phase) / L )
!   where t and phase are in years and L = cycle_length_yr.
!   phase = 0 means the run starts at solar minimum.
!
! Ap is constant (ap_quiet) except during an optional geomagnetic storm
! window [storm_start_day, storm_start_day + storm_length_day).

module space_weather
  use constants
  implicit none

  character(len=16) :: solar_mode = 'constant'
  real(dp) :: f107_const       = 150.0_dp   ! sfu
  real(dp) :: f107_min         = 70.0_dp
  real(dp) :: f107_max         = 220.0_dp
  real(dp) :: cycle_length_yr  = 11.0_dp
  real(dp) :: cycle_phase_yr   = 0.0_dp

  real(dp) :: ap_quiet         = 15.0_dp
  real(dp) :: storm_start_day  = -1.0_dp    ! negative = no storm
  real(dp) :: storm_length_day = 0.0_dp
  real(dp) :: storm_ap         = 200.0_dp

contains

  function f107_at(t_day) result(f107)
    real(dp), intent(in) :: t_day
    real(dp) :: f107, t_yr

    select case (trim(solar_mode))
    case ('constant')
      f107 = f107_const
    case ('cycle')
      t_yr = t_day / day_per_yr
      f107 = f107_min + (f107_max - f107_min) * &
             sin(pi * (t_yr + cycle_phase_yr) / cycle_length_yr)**2
    case default
      print *, 'ERROR: unknown solar_mode: ', trim(solar_mode)
      stop 1
    end select
  end function f107_at

  function ap_at(t_day) result(ap)
    real(dp), intent(in) :: t_day
    real(dp) :: ap

    ap = ap_quiet
    if (storm_start_day >= 0.0_dp) then
      if (t_day >= storm_start_day .and. &
          t_day <  storm_start_day + storm_length_day) ap = storm_ap
    end if
  end function ap_at

end module space_weather
