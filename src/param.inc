
      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 3000, n_bin_max=30, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793238d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      real*8
     $  gmb

      parameter
     $  (gmb = 1.d0+1.d0/6023600.0d0+1.d0/408523.71d0
     $  +1.d0/328900.56d0+1.d0/3098708.0d0+1.d0/1047.3486d0
     $  +1.d0/3497.898d0+1.d0/22902.98d0+1.d0/19412.24d0
     $  +1.d0/1.35d8)
