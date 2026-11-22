! density_profile.f90
! Writes vertical density profiles (150-700 km) for different levels of
! solar and geomagnetic activity, so they can be plotted in Python.
!
! Usage:  density_profile  <output_csv>

program density_profile
  use constants
  use atmosphere
  implicit none

  character(len=256) :: out_file
  integer :: u_out, i
  real(dp) :: h

  if (command_argument_count() /= 1) then
    print *, 'usage: density_profile <output_csv>'
    stop 1
  end if
  call get_command_argument(1, out_file)

  open(newunit=u_out, file=trim(out_file), status='replace', action='write')
  write(u_out, '(a)') 'alt_km,static,f107_70,f107_150,f107_230,f107_150_storm'
  do i = 0, 110
    h = 150.0_dp + 5.0_dp * i
    write(u_out, '(f7.1,5(",",es13.6))') h, static_density(h), &
         sws_density(h,  70.0_dp,  15.0_dp), &
         sws_density(h, 150.0_dp,  15.0_dp), &
         sws_density(h, 230.0_dp,  15.0_dp), &
         sws_density(h, 150.0_dp, 200.0_dp)
  end do
  close(u_out)

  print '(a)', 'Exospheric temperature [K] (Ap = 15):'
  print '(a,f7.1)', '  F10.7 =  70 : ', sws_temperature( 70.0_dp, 15.0_dp)
  print '(a,f7.1)', '  F10.7 = 150 : ', sws_temperature(150.0_dp, 15.0_dp)
  print '(a,f7.1)', '  F10.7 = 230 : ', sws_temperature(230.0_dp, 15.0_dp)
  print '(a,f7.1)', '  storm (F10.7 = 150, Ap = 200) : ', sws_temperature(150.0_dp, 200.0_dp)
end program density_profile
