module datadec

  use poly_dec

  ! define length of array parameters
  integer, parameter :: n_sur_max = 200, n_bin_max=30, n_r_max=10, &
       nw_max = 10

  ! define some useful constants
  real (kind=8), parameter :: Pi = 3.141592653589793238d0, drad = Pi/180.0D0, &
       TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14
  real (kind=8), parameter :: gmb = 1.d0+1.d0/6023600.0d0+1.d0/408523.71d0 &
       +1.d0/328900.56d0+1.d0/3098708.0d0+1.d0/1047.3486d0+1.d0/3497.898d0 &
       +1.d0/22902.98d0+1.d0/19412.24d0+1.d0/1.35d8

  ! define data type to represent survey efficiency and pointings, and objects
  type t_ratecut
     real (kind=8) :: min, max, angle, hwidth
  end type t_ratecut

  type t_orb_m
     real (kind=8) :: a, e, inc, node, peri, m
  end type t_orb_m

  type t_orb_p
     real (kind=8) :: a, e, inc, node, peri, tperi
  end type t_orb_p

  type t_v3d
     real (kind=8) :: x, y, z
  end type t_v3d

  type t_obspos
     real (kind=8) :: jday, r
     type(t_v3d) :: pos
  end type t_obspos

  type t_eff_r
     real (kind=8) :: min, max, mag_lim
     integer :: n
     real (kind=8), dimension(n_bin_max) :: b, e
  end type t_eff_r

  type t_charact
     type(t_ratecut) :: r_cut
     real (kind=8) :: mag_er(6), photf(3), track(3)
     integer :: f, nr
     type(t_eff_r), dimension(n_r_max) :: eff_p
  end type t_charact

  type t_pointing
     real (kind=8) :: ff
     integer :: code
     character :: efnam*80
     type(t_obspos) :: o_pos(2)
     type(t_polygon) :: poly
     type(t_charact) :: c
  end type t_pointing
end module datadec
