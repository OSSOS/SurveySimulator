module datadec

  use parameters
  use poly_dec

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
     character :: f
     integer :: nr
     type(t_eff_r), dimension(n_r_max) :: eff_p
  end type t_charact

  ! Pointing: one FoV row from pointings.list.
  ! Survey: characterization directory basename (unique across loaded surveys).
  ! Block: .eff stem within a survey (unique inside that survey only).
  ! Detection key = trim(survey)//'/'//trim(block); '/' is a key delimiter,
  ! not a filesystem path — open files via eff_file under the survey directory.
  type t_pointing
     real (kind=8) :: ff
     integer :: code
     character(name_len) :: survey
     character(block_len) :: block
     character(eff_name_len) :: eff_file
     character(key_len) :: key
     type(t_obspos) :: o_pos(2)
     type(t_polygon) :: poly
     type(t_charact) :: c
  end type t_pointing

contains 


end module datadec
