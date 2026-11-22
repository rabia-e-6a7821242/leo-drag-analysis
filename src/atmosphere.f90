! atmosphere.f90
! Two simple density models for the thermosphere:
!
!  1) "static"  : piecewise exponential atmosphere from Vallado (2013),
!                 Fundamentals of Astrodynamics and Applications, 4th ed.,
!                 Ch. 8 (exponential atmospheric model table). Values were
!                 cross-checked with the SatelliteToolboxAtmosphericModels.jl
!                 implementation. It does NOT depend on solar activity, so it
!                 is used as a reference.
!
!  2) "sws"     : semi-empirical model from IPS Radio and Space Services /
!                 Bureau of Meteorology note "Satellite Orbital Decay
!                 Calculations", section 2 (equations copied as given there).
!                 The note recommends a 90-day mean F10.7 as input.
!                 Exospheric temperature is driven by solar EUV (F10.7 proxy)
!                 and geomagnetic activity (Ap):
!                     T   = 900 + 2.5 (F10.7 - 70) + 1.5 Ap        [K]
!                     m   = 27 - 0.012 (h - 200)                   [amu]
!                     H   = T / m                                   [km]
!                     rho = 6e-10 exp( -(h - 175) / H )             [kg/m^3]
!                 Valid roughly between 180 and 500 km. Below 180 km we
!                 switch to the static table (the two agree within ~10 %
!                 there, because all SWS profiles pass through 6e-10 at 175 km).

module atmosphere
  use constants
  implicit none
  private
  public :: static_density, sws_density, get_density, sws_temperature

  integer, parameter :: ntab = 28

  ! base altitude [km]
  real(dp), parameter :: h_tab(ntab) = [ &
       0._dp,  25._dp,  30._dp,  40._dp,  50._dp,  60._dp,  70._dp, &
      80._dp,  90._dp, 100._dp, 110._dp, 120._dp, 130._dp, 140._dp, &
     150._dp, 180._dp, 200._dp, 250._dp, 300._dp, 350._dp, 400._dp, &
     450._dp, 500._dp, 600._dp, 700._dp, 800._dp, 900._dp, 1000._dp ]

  ! density at base altitude [kg/m^3]
  real(dp), parameter :: rho_tab(ntab) = [ &
     1.225_dp,    3.899e-2_dp, 1.774e-2_dp, 3.972e-3_dp, 1.057e-3_dp, &
     3.206e-4_dp, 8.770e-5_dp, 1.905e-5_dp, 3.396e-6_dp, 5.297e-7_dp, &
     9.661e-8_dp, 2.438e-8_dp, 8.484e-9_dp, 3.845e-9_dp, 2.070e-9_dp, &
     5.464e-10_dp,2.789e-10_dp,7.248e-11_dp,2.418e-11_dp,9.518e-12_dp, &
     3.725e-12_dp,1.585e-12_dp,6.967e-13_dp,1.454e-13_dp,3.614e-14_dp, &
     1.170e-14_dp,5.245e-15_dp,3.019e-15_dp ]

  ! scale height [km]
  real(dp), parameter :: sh_tab(ntab) = [ &
     7.249_dp,  6.349_dp,  6.682_dp,  7.554_dp,  8.382_dp,  7.714_dp, &
     6.549_dp,  5.799_dp,  5.382_dp,  5.877_dp,  7.263_dp,  9.473_dp, &
    12.636_dp, 16.149_dp, 22.523_dp, 29.740_dp, 37.105_dp, 45.546_dp, &
    53.628_dp, 53.298_dp, 58.515_dp, 60.828_dp, 63.822_dp, 71.835_dp, &
    88.667_dp,124.64_dp, 181.05_dp, 268.00_dp ]

contains

  ! ------------------------------------------------------------------
  function static_density(h_km) result(rho)
    real(dp), intent(in) :: h_km
    real(dp) :: rho
    integer :: i

    i = 1
    do while (i < ntab)
      if (h_km < h_tab(i+1)) exit
      i = i + 1
    end do
    rho = rho_tab(i) * exp(-(h_km - h_tab(i)) / sh_tab(i))
  end function static_density

  ! ------------------------------------------------------------------
  function sws_temperature(f107, ap) result(t_exo)
    real(dp), intent(in) :: f107, ap
    real(dp) :: t_exo
    t_exo = 900.0_dp + 2.5_dp * (f107 - 70.0_dp) + 1.5_dp * ap
  end function sws_temperature

  ! ------------------------------------------------------------------
  function sws_density(h_km, f107, ap) result(rho)
    real(dp), intent(in) :: h_km, f107, ap
    real(dp) :: rho, t_exo, m_mol, h_scale

    if (h_km < 180.0_dp) then
      rho = static_density(h_km)
      return
    end if

    t_exo   = sws_temperature(f107, ap)
    m_mol   = max(27.0_dp - 0.012_dp * (h_km - 200.0_dp), 16.0_dp)  ! never lighter than atomic O
    h_scale = t_exo / m_mol
    rho     = 6.0e-10_dp * exp(-(h_km - 175.0_dp) / h_scale)
  end function sws_density

  ! ------------------------------------------------------------------
  function get_density(h_km, f107, ap, model) result(rho)
    real(dp), intent(in) :: h_km, f107, ap
    character(len=*), intent(in) :: model
    real(dp) :: rho

    select case (trim(model))
    case ('static')
      rho = static_density(h_km)
    case ('sws')
      rho = sws_density(h_km, f107, ap)
    case default
      print *, 'ERROR: unknown density model: ', trim(model)
      stop 1
    end select
  end function get_density

end module atmosphere
