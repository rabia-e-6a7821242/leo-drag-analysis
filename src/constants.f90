! constants.f90
! Physical constants used in the whole project (SI units unless noted).

module constants
  use iso_fortran_env, only: dp => real64
  implicit none

  real(dp), parameter :: pi          = 3.14159265358979323846_dp
  real(dp), parameter :: deg2rad     = pi / 180.0_dp
  real(dp), parameter :: mu_earth    = 3.986004418e14_dp   ! m^3/s^2, GM of Earth
  real(dp), parameter :: r_earth     = 6378.137e3_dp       ! m, WGS-84 equatorial radius
  real(dp), parameter :: sec_per_day = 86400.0_dp
  real(dp), parameter :: day_per_yr  = 365.25_dp

end module constants
