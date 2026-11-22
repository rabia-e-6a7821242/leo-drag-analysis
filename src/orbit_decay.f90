! orbit_decay.f90
! Orbital decay of a LEO satellite on a (near-)circular orbit due to drag.
!
! For a circular orbit the drag-induced change of semi-major axis is
!     da/dt = - rho * (Cd * A / m) * sqrt(mu * a)
! and the satellite moves along the orbit with the mean motion
!     du/dt = sqrt(mu / a^3)
! The system (a, u) is integrated with classical 4th-order Runge-Kutta.
! u is needed only because the frontal area of the solar arrays depends
! on the position along the orbit (see spacecraft.f90).
!
! Usage:  orbit_decay  <namelist_file>  <output_csv>

program orbit_decay
  use constants
  use atmosphere,    only: get_density
  use space_weather
  use spacecraft
  implicit none

  ! ---- run settings (read from namelist) ----
  real(dp) :: h0_km = 450.0_dp, h_end_km = 150.0_dp
  character(len=16) :: density_model = 'sws'
  real(dp) :: dt_s = 60.0_dp, out_dt_s = 21600.0_dp, t_max_day = 5475.0_dp

  namelist /satellite/   mass_kg, cd, area_body, area_panel, area_mode, &
                         beta_amp_deg, beta_period_day
  namelist /orbit/       h0_km, h_end_km
  namelist /environment/ density_model, solar_mode, f107_const, f107_min, &
                         f107_max, cycle_length_yr, cycle_phase_yr, ap_quiet, &
                         storm_start_day, storm_length_day, storm_ap
  namelist /run/         dt_s, out_dt_s, t_max_day

  character(len=256) :: nml_file, out_file
  integer :: u_nml, u_out, ios
  real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2)
  real(dp) :: t, t_next_out, t_day, h_km
  real(dp) :: sum_area, sum_rho, rho_now, area_now
  integer  :: n_avg
  integer(8) :: nstep
  logical :: reentered

  if (command_argument_count() /= 2) then
    print *, 'usage: orbit_decay <namelist_file> <output_csv>'
    stop 1
  end if
  call get_command_argument(1, nml_file)
  call get_command_argument(2, out_file)

  ! ---- read input; each group is optional ----
  open(newunit=u_nml, file=trim(nml_file), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'ERROR: cannot open ', trim(nml_file)
    stop 1
  end if
  read(u_nml, nml=satellite,   iostat=ios); rewind(u_nml)
  read(u_nml, nml=orbit,       iostat=ios); rewind(u_nml)
  read(u_nml, nml=environment, iostat=ios); rewind(u_nml)
  read(u_nml, nml=run,         iostat=ios)
  close(u_nml)

  ! ---- initial state ----
  y(1) = r_earth + h0_km * 1.0e3_dp   ! semi-major axis [m]
  y(2) = 0.0_dp                       ! argument from Sun direction [rad]
  t = 0.0_dp
  t_next_out = 0.0_dp
  sum_area = 0.0_dp; sum_rho = 0.0_dp; n_avg = 0
  nstep = 0
  reentered = .false.

  open(newunit=u_out, file=trim(out_file), status='replace', action='write')
  write(u_out, '(a)') 'time_day,alt_km,f107,ap,density_kgm3,area_m2,beta_deg,period_min'

  do
    t_day = t / sec_per_day
    h_km  = (y(1) - r_earth) / 1.0e3_dp

    ! values at this step, averaged over the output interval
    rho_now  = get_density(h_km, f107_at(t_day), ap_at(t_day), density_model)
    area_now = frontal_area(y(2), t_day)
    sum_rho  = sum_rho + rho_now
    sum_area = sum_area + area_now
    n_avg    = n_avg + 1

    if (h_km <= h_end_km) reentered = .true.

    if (t >= t_next_out .or. reentered .or. t_day >= t_max_day) then
      call write_row()
      t_next_out = t_next_out + out_dt_s
    end if

    if (reentered .or. t_day >= t_max_day) exit

    ! ---- RK4 step ----
    k1 = deriv(t,                y)
    k2 = deriv(t + 0.5_dp*dt_s,  y + 0.5_dp*dt_s*k1)
    k3 = deriv(t + 0.5_dp*dt_s,  y + 0.5_dp*dt_s*k2)
    k4 = deriv(t + dt_s,         y + dt_s*k3)
    y  = y + dt_s / 6.0_dp * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
    y(2) = modulo(y(2), 2.0_dp * pi)
    t  = t + dt_s
    nstep = nstep + 1
  end do
  close(u_out)

  ! summary for the Python driver (parsed line by line)
  write(*, '(a,a)')      'CASE          = ', trim(nml_file)
  write(*, '(a,l1)')     'REENTERED     = ', reentered
  write(*, '(a,f12.3)')  'LIFETIME_DAYS = ', t / sec_per_day
  write(*, '(a,f12.3)')  'FINAL_ALT_KM  = ', (y(1) - r_earth) / 1.0e3_dp
  write(*, '(a,i0)')     'STEPS         = ', nstep

contains

  function deriv(time, state) result(dydt)
    real(dp), intent(in) :: time, state(2)
    real(dp) :: dydt(2), a, alt, rho, area, td

    a    = state(1)
    alt  = (a - r_earth) / 1.0e3_dp
    td   = time / sec_per_day
    rho  = get_density(alt, f107_at(td), ap_at(td), density_model)
    area = frontal_area(state(2), td)

    dydt(1) = -rho * (cd * area / mass_kg) * sqrt(mu_earth * a)
    dydt(2) = sqrt(mu_earth / a**3)
  end function deriv

  subroutine write_row()
    real(dp) :: period_min
    period_min = 2.0_dp * pi * sqrt(y(1)**3 / mu_earth) / 60.0_dp
    write(u_out, '(f12.5,",",f10.4,",",f8.2,",",f8.2,",",es13.6,",",f8.5,",",f8.3,",",f9.4)') &
         t_day, h_km, f107_at(t_day), ap_at(t_day), sum_rho / n_avg, &
         sum_area / n_avg, beta_at(t_day) / deg2rad, period_min
    sum_rho = 0.0_dp; sum_area = 0.0_dp; n_avg = 0
  end subroutine write_row

end program orbit_decay
