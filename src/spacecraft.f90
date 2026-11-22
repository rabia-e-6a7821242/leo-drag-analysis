! spacecraft.f90
! Satellite properties and the frontal (ram-facing) area model.
!
! The satellite is a box body plus sun-tracking solar arrays.
! The body always shows area_body to the flow. The arrays point their
! normal at the Sun, so the part of the array area that faces the
! velocity vector changes along the orbit.
!
! Geometry (circular orbit, Sun at angle beta above the orbit plane):
!   s_hat = cos(beta) x_hat + sin(beta) z_hat         (z = orbit normal)
!   v_hat = -sin(u) x_hat + cos(u) y_hat              (u = angle from Sun dir.)
!   projected array area = area_panel * |s_hat . v_hat|
!                        = area_panel * |cos(beta) sin(u)|
!
! area_mode:
!   'fixed_max' : body + full panel area all the time (worst case)
!   'fixed_avg' : body + orbit-averaged panel area for beta = 0, (2/pi)*panel
!   'varying'   : instantaneous area from the formula above

module spacecraft
  use constants
  implicit none

  real(dp) :: mass_kg         = 150.0_dp
  real(dp) :: cd              = 2.2_dp
  real(dp) :: area_body       = 0.6_dp    ! m^2
  real(dp) :: area_panel      = 2.4_dp    ! m^2 (both wings)
  character(len=16) :: area_mode = 'varying'

  ! The beta angle slowly changes because of orbit precession and the
  ! Earth's motion around the Sun. Simple sinusoidal model:
  real(dp) :: beta_amp_deg    = 0.0_dp
  real(dp) :: beta_period_day = 60.0_dp

contains

  function beta_at(t_day) result(beta)
    real(dp), intent(in) :: t_day
    real(dp) :: beta
    beta = beta_amp_deg * deg2rad * sin(2.0_dp * pi * t_day / beta_period_day)
  end function beta_at

  function frontal_area(u, t_day) result(area)
    real(dp), intent(in) :: u, t_day
    real(dp) :: area

    select case (trim(area_mode))
    case ('fixed_max')
      area = area_body + area_panel
    case ('fixed_avg')
      area = area_body + (2.0_dp / pi) * area_panel
    case ('varying')
      area = area_body + area_panel * abs(cos(beta_at(t_day)) * sin(u))
    case default
      print *, 'ERROR: unknown area_mode: ', trim(area_mode)
      stop 1
    end select
  end function frontal_area

end module spacecraft
