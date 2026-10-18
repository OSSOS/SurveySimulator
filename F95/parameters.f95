module parameters

  ! define length of array parameters
  integer, parameter :: n_sur_max = 2000, n_bin_max=30, n_r_max=10, &
       nw_max = 10

  ! Filesystem path and survey/block identity lengths.
  ! path_len: sole size for any filesystem path passed into the code.
  ! name_len / block_len: halves of the detection key (survey and block).
  ! key_len: composite "survey/block" written to detect output ( '/' is a
  ! key delimiter, not a path separator — do not open key as a file).
  integer, parameter :: path_len  = 2048
  integer, parameter :: name_len  = 32
  integer, parameter :: block_len = 32
  integer, parameter :: key_len   = name_len + 1 + block_len
  ! Efficiency basename as listed in pointings.list (stem + optional ".eff").
  integer, parameter :: eff_name_len = block_len + 4

  ! define some useful constants
  real (kind=8), parameter :: Pi = 3.141592653589793238d0, drad = Pi/180.0D0, &
       TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14
  real (kind=8), parameter :: gmb = 1.d0+1.d0/6023600.0d0+1.d0/408523.71d0 &
       +1.d0/328900.56d0+1.d0/3098708.0d0+1.d0/1047.3486d0+1.d0/3497.898d0 &
       +1.d0/22902.98d0+1.d0/19412.24d0+1.d0/1.35d8
  ! define the order of filters in classic fortran model in/out files
  character(len=10), parameter :: filters = "grizuVBRIw"
  integer, parameter :: number_of_predefined_filters = 10

end module parameters
