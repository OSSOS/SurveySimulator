module xvutils

  use elemutils
  use rot
  use ioutils

contains
  subroutine BaryXV (jday, pos, vel, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives the position and velocity in ecliptic heliocentric
! reference frame of the Solar System barycenter at a given time. Uses
! PlanetXV.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     jday  : Time of elements (Julian day) (R8)
!
! OUTPUT
!     pos   : Position vector (3*R8)
!     vel   : Velocity vector (3*R8)
!     ierr  : Error code (I4)
!                0 : nominal run
!               20 : jday out of range
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) jday
!f2py intent(out) pos
!f2py intent(out) vel
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(out) :: pos, vel
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: jday
    integer, parameter :: n_planets = 9
    real (kind=8), parameter :: jday_min = 2415020.0, jday_max = 2488070.0, &
         mu = TwoPi**2
    type(t_v3d) :: pos_p, vel_p
    integer :: ind, istat
    real (kind=8) :: masses(n_planets), msys, mp

    data masses /6023600.0d0,   &  ! Mercury
                  408523.71d0,  &  ! Venus
                  328900.56d0,  &  ! Earth + Moon
                 3098708.0d0,   &  ! Mars
                    1047.3486d0,&  ! Jupiter
                    3497.898d0, &  ! Saturn
                   22902.98d0,  &  ! Uranus
                   19412.24d0,  &  ! Neptune
                       1.35d8/     ! Pluto

    ierr = 0

! Jday out of range.
    if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
       ierr = 20
       return
    end if

! Ok, do the math.
    pos%x = 0.d0
    pos%y = 0.d0
    pos%z = 0.d0
    vel%x = 0.d0
    vel%y = 0.d0
    vel%z = 0.d0
    msys = mu
    do ind = 1, n_planets
       call PlanetXV (ind, jday, pos_p, vel_p, istat)
       if (istat .ne. 0) then
          ierr = istat
          return
       end if
       mp = mu/masses(ind)
       msys = msys + mp
       pos%x = pos%x + mp*pos_p%x
       pos%y = pos%y + mp*pos_p%y
       pos%z = pos%z + mp*pos_p%z
       vel%x = vel%x + mp*vel_p%x
       vel%y = vel%y + mp*vel_p%y
       vel%z = vel%z + mp*vel_p%z
    end do
    pos%x = pos%x/msys
    pos%y = pos%y/msys
    pos%z = pos%z/msys
    vel%x = vel%x/msys
    vel%y = vel%y/msys
    vel%z = vel%z/msys

    return

  end subroutine BaryXV

  subroutine DistSunEcl (jday, pos, r)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes distance from object defined by barycentric
! ecliptic coordinates to Sun at given jday.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     jday  : Time (R8)
!     pos   : Object barycentric ecliptic cartsian coordinates (3*R8)
!
! OUTPUT
!     r     : Distance from object to Sun (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) jday
!f2py intent(in) pos
!f2py intent(out) r
    implicit none

    type(t_v3d), intent(in) :: pos
    real (kind=8), intent(in) :: jday
    real (kind=8), intent(out) :: r
    type(t_v3d), save :: pos_b, vel_b
    integer :: ierr
    real (kind=8), save :: jday_b
    logical, save :: done

    data done /.false./

    if ((.not. done) .or. (jday .ne. jday_b)) then
       call BaryXV (jday, pos_b, vel_b, ierr)
       if (ierr .ne. 0) then
          write (6, *) 'Problem while getting barycenter position'
       end if
       done = .true.
       jday_b = jday
    end if

    r = dsqrt((pos%x + pos_b%x)**2 + (pos%y + pos_b%y)**2 &
         + (pos%z + pos_b%z)**2)
    return

  end subroutine DistSunEcl

  subroutine ObsPos (code, t, pos, vel, r, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the cartesian coordinates of the observatory.
! Reference frame : ICRF
! Units : AU
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2003
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     code  : Observatory code (I4)
!              001 : GAIA
!              002 : Geocentric, Mignard's code
!              500 : Geocentric
!              <0  : -ve of a file handle to read JPL state vectors from
!     t     : Time of observation (Julian day, not MJD) (R8)
!
! OUTPUT
!     pos   : Cartesian coordinates of observatory (AU) *(R8)
!     vel   : Cartesian velocity of observatory (AU) *(R8)
!     r     : Distance from observatory to Sun (AU) (R8)
!     ierr  : Error code (I4)
!                0 : nominal run
!               10 : unknown observatory code
!              100 : date of call earlier than xjdbeg
!              200 : date of call later   than xjdend
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) code
!f2py intent(in) t
!f2py intent(out) pos
!f2py intent(out) vel
!f2py intent(out) r
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(out) :: pos, vel
    integer, intent(in) :: code
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: t
    real (kind=8), intent(out) :: r
    type(t_v3d) :: pos_b, vel_b
    integer :: istat
    real (kind=8), parameter :: km2AU = 149597870.691d0
    real (kind=8) :: v_pos(3)

    ierr = 0
    if (code .eq. 1) then
       ierr = 10
       return
    else if (code .eq. 2) then
       ierr = 10
       return 
    else if (code .lt. 0) then
! This is  a file handle, send to read_jpl_csv
       call read_jpl_csv(-1*code, t, pos, vel, ierr)
       if (ierr .ne. 0) then
          write(0,*) "Failed while reading JPL Ephem for time ",t," using LUN: ",-1*code
          return 
       end if
    else
       
! Get heliocentric position of Earth.
       call newcomb (t, v_pos)
       pos%x = -v_pos(1)
       pos%y = -v_pos(2)
       pos%z = -v_pos(3)
       vel%x = 0.d0
       vel%y = 0.d0
       vel%z = 0.d0

! Now move to barycenter. First get barycenter position (ecliptic).
       call BaryXV (t, pos_b, vel_b, istat)
       if (ierr .ne. 0) then
          write (6, *) 'Problem while getting barycenter position'
       end if

! Convert barycenter position to Equatorial.
       call equat_ecl (-1, pos_b, pos_b, ierr)
       if (ierr .ne. 0) then
          write (6, *) 'Problem in conversion ecliptic -> equatorial'
       end if
       call equat_ecl (-1, vel_b, vel_b, ierr)
       if (ierr .ne. 0) then
          write (6, *) 'Problem in conversion ecliptic -> equatorial'
       end if
       pos%x = pos%x - pos_b%x
       pos%y = pos%y - pos_b%y
       pos%z = pos%z - pos_b%z
       vel%x = vel%x - vel_b%x
       vel%y = vel%y - vel_b%y
       vel%z = vel%z - vel_b%z
       if (code .eq. 500) then
       else
          ierr = 10
          return
       end if
    end if
! Finally, computes distance from observatory to Sun.
    r = dsqrt((pos%x + pos_b%x)**2 + (pos%y + pos_b%y)**2 &
         + (pos%z + pos_b%z)**2)

  end subroutine ObsPos

  subroutine PlanetElem (ind, jday, o_m, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives the osculating elements in heliocentric reference
! frame of a planet at a given time. From given elements and rates.
! Valid roughly from 1900 to 2100.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : September 2019
!             Added Pluto (index 9) and used new orbital elements from
!             JPL.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     ind   : Planet index (I4)
!                1 : Mercury
!                2 : Venus
!                3 : Earth
!                4 : Mars
!                5 : Jupiter
!                6 : Saturn
!                7 : Uranus
!                8 : Neptune
!                9 : Pluto
!     jday  : Time of elements (Julian day) (R8)
!
! OUTPUT
!     o_m   : Orbital elements (t_orb_m)
!     ierr  : Error code (I4)
!                0 : nominal run
!               10 : Unknown planet index
!               20 : jday out of range
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ind
!f2py intent(in) jday
!f2py intent(out) o_m
!f2py intent(out) ierr
    implicit none

    type(t_orb_m), intent(out) :: o_m
    integer, intent(in) :: ind
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: jday
    real (kind=8), parameter :: jday_min = 2378496.5d0, jday_max = 2470171.5d0
    integer, parameter :: n_planets = 9
    real (kind=8) :: a_p(n_planets), e_p(n_planets), i_p(n_planets), &
         node_p(n_planets), peri_p(n_planets), capl_p(n_planets), &
         da_p(n_planets), de_p(n_planets), di_p(n_planets), &
         dnode_p(n_planets), dperi_p(n_planets), dcapl_p(n_planets), &
         jday_p, dt

! Elements are given in heliocentric reference frame.
    data &
! Reference time is time of elements, 2451545.0 = 2000-01-01.5
         jday_p /2451545.0d0/, &
! Semi-major axis (au)    
         a_p / 0.38709927d0, 0.72333566d0, 1.00000261d0, 1.52371034d0, &
         5.20288700d0, 9.53667594d0, 19.18916464d0, 30.06992276d0, &
         39.48211675d0/, &
! Eccentricity
         e_p / 0.20563593d0, 0.00677672d0, 0.01671123d0, 0.09339410d0, &
         0.04838624d0, 0.05386179d0, 0.04725744d0, 0.00859048d0, &
         0.24882730d0/, &
! Ecliptic inclination (deg)
         i_p / 7.00497902d0, 3.39467605d0, -0.00001531d0, 1.84969142d0, &
         1.30439695d0, 2.48599187d0, 0.77263783d0, 1.77004347d0, &
         17.14001206d0/, &
! Longitude of node (deg)
         node_p / 48.33076593d0, 76.67984255d0, 0.0d0, 49.55953891d0, &
         100.47390909d0, 113.66242448d0, 74.01692503d0, 131.78422574d0, &
         110.30393684d0/, &
! Longitude of pericenter (deg)
         peri_p / 77.45779628d0, 131.60246718d0, 102.93768193d0, &
         -23.94362959d0, 14.72847983d0, 92.59887831d0, 170.95427630d0, &
         44.96476227d0, 224.06891629d0/, &
! Mean longitude (deg)
         capl_p / 252.25032350d0, 181.97909950d0, 100.46457166d0, &
         -4.55343205d0, 34.39644051d0, 49.95424423d0, 313.23810451d0, &
         -55.12002969d0, 238.92903833d0/, &
! Rate is in au per century
         da_p / 0.00000037d0, 0.00000390d0, 0.00000562d0, 0.00001847d0, &
         -0.00011607d0, -0.00125060d0, -0.00196176d0, 0.00026291d0, &
         -0.00031596d0/, &
! Rate is in 1 per century
         de_p / 0.00001906d0, -0.00004107d0, -0.00004392d0, 0.00007882d0, &
         -0.00013253d0, -0.00050991d0, -0.00004397d0, 0.00005105d0, &
         0.00005170d0/, &
! Rates below are given in degrees per century
         di_p / -0.00594749d0, -0.00078890d0, -0.01294668d0, &
         -0.00813131d0, -0.00183714d0, 0.00193609d0, -0.00242939d0, &
         0.00035372d0, 0.00004818d0/, &
         dnode_p / -0.12534081d0, -0.27769418d0, 0.0d0, -0.29257343d0, &
         0.20469106d0, -0.28867794d0, 0.04240589d0, -0.00508664d0, &
         -0.01183482d0/, &
         dperi_p / 0.16047689d0, 0.00268329d0, 0.32327364d0, &
         0.44441088d0, 0.21252668d0, -0.41897216d0, 0.40805281d0, &
         -0.32241464d0, -0.04062942d0/, &
         dcapl_p / 149472.67411175d0, 58517.81538729d0, 35999.37244981d0, &
         19140.30268499d0, 3034.74612775d0, 1222.49362201d0, &
         428.48202785d0, 218.45945325d0, 145.20780515d0/

    ierr = 0

! Wrong planet index.
    if ((ind .lt. 1) .or. (ind .gt. n_planets)) then
       ierr = 10
       return
    end if

! Jday out of range.
    if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
       ierr = 20
       return
    end if

! Ok, do the math.
! Rates are given in degrees per century, and angles in degrees. Also we
! have all longitudes, when we need arguments. So tranform in arguments
! and radians.
    dt = (jday - jday_p)/36525.d0
    o_m%a = a_p(ind) + dt*da_p(ind)
    o_m%e = e_p(ind) + dt*de_p(ind)
    o_m%inc = (i_p(ind) + dt*di_p(ind))*drad
    o_m%node = (node_p(ind) + dt*dnode_p(ind))*drad

! Longitude of pericenter
    o_m%peri = (peri_p(ind) + dt*dperi_p(ind))*drad

! Mean longitude. Change to mean anomaly
    o_m%m = (capl_p(ind) + dt*dcapl_p(ind))*drad - o_m%peri

! Now get argument of pericenter
    o_m%peri = o_m%peri - o_m%node

    return

  end subroutine PlanetElem

  subroutine PlanetXV (ind, jday, pos, vel, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives the position and velocity in ecliptic heliocentric
! reference frame of a planet at a given time. Uses PlanetElem.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : January 2019    
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     ind   : Planet index (I4)
!                1 : Mercury
!                2 : Venus
!                3 : Earth
!                4 : Mars
!                5 : Jupiter
!                6 : Saturn
!                7 : Uranus
!                8 : Neptune
!                9 : Pluto    
!     jday  : Time of elements (Julian day) (R8)
!
! OUTPUT
!     pos   : Position vector (3*R8)
!     vel   : Velocity vector (3*R8)
!     ierr  : Error code (I4)
!                0 : nominal run
!               10 : Unknown planet index
!               20 : jday out of range
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ind
!f2py intent(in) jday
!f2py intent(out) pos
!f2py intent(out) vel
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(out) :: pos, vel
    integer, intent(in) :: ind
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: jday
    real (kind=8), parameter :: jday_min = 2415020.0d0, jday_max = 2488070.0d0, &
         mu = TwoPi**2
    integer, parameter :: n_planets = 9
    integer :: istat
    type(t_orb_m) :: o_m
    real (kind=8) :: masses(n_planets), gm

    data masses /6023600.0d0,   &  ! Mercury
                  408523.7d0,   &  ! Venus
                  328900.56d0,  &  ! Earth + Moon
                 3098708.0d0,   &  ! Mars
                    1047.3486d0,&  ! Jupiter
                    3497.898d0, &  ! Saturn
                   22902.98d0,  &  ! Uranus
                   19412.24d0,  &  ! Neptune
                       1.35d8/     ! Pluto

    ierr = 0

! Wrong planet index.
    if ((ind .lt. 1) .or. (ind .gt. n_planets)) then
       ierr = 10
       return
    end if

! Jday out of range.
    if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
       ierr = 20
       return
    end if

! Ok, do the math.
    call PlanetElem (ind, jday, o_m, istat)
    if (istat .ne. 0) then
       ierr = istat
       return
    end if
    gm = mu*(1.d0 + 1.d0/masses(ind))
    call coord_cart (gm, o_m, pos, vel)

    return

  end subroutine PlanetXV

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives the cartesian coordinates of the Earth at a given
! time in a heliocentric, equatorial coordinate system. This is a self
! -contained subroutine. One can also use a progamme that reads the
! JPL ephemerides DE405 (see ss_state.f program).
!
! Tested at different times, clearly gives equatorial coordinates.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  subroutine newcomb (DJ,ROUT)
!f2py intent(in) dj
!f2py intent(out) rout
    implicit none

    real (kind=8), intent(in) :: dj
    real (kind=8), intent(out) :: rout(3)
    real (kind=8) :: Z(8,121),Z0(80),Z1(80),Z2(80),Z3(80),Z4(80),Z5(80),Z6(80), &
         Z7(80),Z8(80),Z9(80),ZA(80),ZB(80),ZC(8),ROU(3),CC(3,3), &
         EL(21),C(5),RIN(3),X(3),Y(3),W(3), TWOPI, STR, RTD, OBL0, &
         C1, C2, C3, TFIN, T, TP, T18, V, V2, V3, T50, D, ARG, DL, DR, DBLARG, &
         DB, GME, GV, GMA, GJ, GS, RO, CS, SS, BS, BC, BR, RS, RC, OBL, T2, T3, &
         RDPSC, DPTTY, SDL, DT, DT2, DY, DY2, DY3, AK, AW, BN, SINK, COSK, &
         SINN, COSN, RIH, RIG, DIG, SINW, COSW, DG
    integer :: i, k

    EQUIVALENCE(Z(1,  1),Z0(1))
    EQUIVALENCE(Z(1, 11),Z1(1))
    EQUIVALENCE(Z(1, 21),Z2(1))
    EQUIVALENCE(Z(1, 31),Z3(1))
    EQUIVALENCE(Z(1, 41),Z4(1))
    EQUIVALENCE(Z(1, 51),Z5(1))
    EQUIVALENCE(Z(1, 61),Z6(1))
    EQUIVALENCE(Z(1, 71),Z7(1))
    EQUIVALENCE(Z(1, 81),Z8(1))
    EQUIVALENCE(Z(1, 91),Z9(1))
    EQUIVALENCE(Z(1,101),ZA(1))
    EQUIVALENCE(Z(1,111),ZB(1))
    EQUIVALENCE(Z(1,121),ZC(1))
    DATA TWOPI/6.283185307179586D0/
    DATA STR/206264806.2470964D0/
    DATA RTD/57.29577951308232D0/
    DATA OBL0,C1,C2,C3/23.4522944D0,.1301250D-1,.163889D-5,.50278D-6/
    DATA TFIN/2433282.5D0/
    DATA Z0/   - 1.,  1., -   6.,-  11.,    26.,-  12.,     0.,    0., &
         - 1.,  2., -   3.,-   3., -   4.,    5.,     0.,    0., &
         - 1.,  3.,    15.,-   1., -   1.,-  18.,     0.,    0., &
         - 1.,  4.,    19.,-  13., -   2.,-   4.,     0.,    0., &
         - 1.,  0.,    33.,-  67., -  85.,-  39.,    24.,-  17., &
         - 1.,  1.,  2350.,-4223., -2059.,-1145., -   4.,    3., &
         - 1.,  2., -  65.,-  34.,    68.,-  14.,     6.,-  92., &
         - 1.,  3., -   3.,-   8.,    14.,-   8.,     1.,    7., &
         - 2.,  0., -   3.,    1.,     0.,    4.,     0.,    0., &
         - 2.,  1., -  99.,   60.,    84.,  136.,    23.,-   3./
    DATA Z1/   - 2.,  2., -4696., 2899.,  3588., 5815.,    10.,-   6., &
           - 2.,  3.,  1793.,-1735., - 595.,- 631.,    37.,-  56., &
           - 2.,  4.,    30.,-  33.,    40.,   33.,     5.,-  13., &
           - 3.,  2., -  13.,    1.,     0.,   21.,    13.,    5., &
           - 3.,  3., - 665.,   27.,    44., 1043.,     8.,    1., &
           - 3.,  4.,  1506.,- 396., - 381.,-1446.,   185.,- 100., &
           - 3.,  5.,   762.,- 683.,   126.,  148.,     6.,-   3., &
           - 3.,  6.,    12.,-  12.,    14.,   13., -   2.,    4., &
           - 4.,  3., -   3.,-   1.,     0.,    6.,     4.,    5., &
           - 4.,  4., - 188.,-  93., - 166.,  337.,     0.,    0./
    DATA Z2/   - 4.,  5., - 139.,-  38., -  51.,  189., -  31.,-   1., &
         - 4.,  6.,   146.,-  42., -  25.,-  91.,    12.,    0., &
         - 4.,  7.,     5.,-   4.,     3.,    5.,     0.,    0., &
         - 5.,  5., -  47.,-  69., - 134.,   93.,     0.,    0., &
         - 5.,  6., -  28.,-  25., -  39.,   43., -   8.,-   4., &
         - 5.,  7., - 119.,-  33., -  37.,  136., -  18.,-   6., &
         - 5.,  8.,   154.,-   1.,     0.,-  26.,     0.,    0., &
         - 6.,  5.,     0.,    0.,     0.,    0., -   2.,    6., &
         - 6.,  6., -   4.,-  38., -  80.,    8.,     0.,    0., &
         - 6.,  7., -   4.,-  13., -  24.,    7., -   2.,-   3./
    DATA Z3/   - 6.,  8., -   6.,-   7., -  10.,   10., -   2.,-   3., &
         - 6.,  9.,    14.,    3.,     3.,-  12.,     0.,    0., &
         - 7.,  7.,     8.,-  18., -  38.,-  17.,     0.,    0., &
         - 7.,  8.,     1.,-   6., -  12.,-   3.,     0.,    0., &
         - 7.,  9.,     1.,-   3., -   4.,    1.,     0.,    0., &
         - 7., 10.,     0.,    0., -   3.,    3.,     0.,    0., &
         - 8.,  8.,     9.,-   7., -  14.,-  19.,     0.,    0., &
         - 8.,  9.,     0.,    0., -   5.,-   4.,     0.,    0., &
         - 8., 12., -   8.,-  41., -  43.,    8., -   5.,-   9., &
         - 8., 13.,     0.,    0., -   9.,-   8.,     0.,    0./
    DATA Z4/   - 8., 14.,    21.,   24., -  25.,   22.,     0.,    0., &
         - 9.,  9.,     6.,-   1., -   2.,-  13.,     0.,    0., &
         - 9., 10.,     0.,    0., -   1.,-   4.,     0.,    0., &
         -10., 10.,     3.,    1.,     3.,-   7.,     0.,    0., &
         1.,- 2., -   5.,-   4., -   5.,    6.,     0.,    0., &
         1.,- 1., - 216.,- 167., -  92.,  119.,     0.,    0., &
         1.,  0., -   8.,-  47.,    27.,-   6.,     0.,    0., &
         2.,- 3.,    40.,-  10., -  13.,-  50.,     0.,    0., &
         2.,- 2.,  1960.,- 566., - 572.,-1973.,     0.,-   8., &
         2.,- 1., -1656.,- 616.,    64.,- 137.,     0.,    0./
    DATA Z5/     2.,  0., -  24.,   15., -  18.,-  25., -   8.,    2., &
         3.,- 4.,     1.,-   4., -   6.,    0.,     0.,    0., &
         3.,- 3.,    53.,- 118., - 154.,-  67.,     0.,    0., &
         3.,- 2.,   395.,- 153., -  77.,- 201.,     0.,    0., &
         3.,- 1.,     8.,    1.,     0.,    6.,     0.,    0., &
         4.,- 4.,    11.,   32.,    46.,-  17.,     0.,    0., &
         4.,- 3., - 131.,  482.,   460.,  125.,     7.,    1., &
         4.,- 2.,   525.,- 256.,    43.,   96.,     0.,    0., &
         4.,- 1.,     7.,-   5.,     6.,    8.,     0.,    0., &
         5.,- 5., -   7.,    1.,     0.,   12.,     0.,    0./
    DATA Z6/     5.,- 4.,    49.,   69.,    87.,-  62.,     0.,    0., &
         5.,- 3., -  38.,  200.,    87.,   17.,     0.,    0., &
         5.,- 2.,     3.,    1., -   1.,    3.,     0.,    0., &
         6.,- 6.,     0.,    0., -   4.,-   3.,     0.,    0., &
         6.,- 5., -  20.,-   2., -   3.,   30.,     0.,    0., &
         6.,- 4., - 104.,- 113., - 102.,   94.,     0.,    0., &
         6.,- 3., -  11.,  100., -  27.,-   4.,     0.,    0., &
         7.,- 6.,     3.,-   5., -   9.,-   5.,     0.,    0., &
         7.,- 5., -  49.,    3.,     4.,   60.,     0.,    0., &
         7.,- 4., -  78.,-  72., -  26.,   28.,     0.,    0./
    DATA Z7/     8.,- 7.,     1.,    3.,     5.,-   1.,     0.,    0., &
         8.,- 6.,     6.,-   8., -  12.,-   9.,     0.,    0., &
         8.,- 5.,    51.,-  10., -   8.,-  44.,     0.,    0., &
         8.,- 4., -  17.,-  12.,     5.,-   6.,     0.,    0., &
         9.,- 7.,     2.,    3.,     5.,-   3.,     0.,    0., &
         9.,- 6.,    13.,-  25., -  30.,-  16.,     0.,    0., &
         9.,- 5.,    60.,-  15., -   4.,-  17.,     0.,    0., &
         10.,- 7.,     2.,    5.,     7.,-   3.,     0.,    0., &
         10.,- 6., -   7.,   18.,    14.,    6.,     0.,    0., &
         10.,- 5.,     5.,-   2.,     0.,    0.,     0.,    0./
    DATA Z8/    11.,- 7.,     9.,   15.,    17.,-  10.,     0.,    0., &
         11.,- 6., -  12.,   42.,     8.,    3.,     0.,    0., &
         12.,- 7., -   4.,-   5., -   4.,    3.,     0.,    0., &
         13.,- 8., -  13.,-   1., -   1.,   15.,     0.,    0., &
         13.,- 7., -  30.,-  33., -   4.,    3.,     0.,    0., &
         15.,- 9.,    13.,-  16., -  17.,-  14.,     0.,    0., &
         15.,- 8.,   200.,-  30., -   1.,-   6.,     0.,    0., &
         17.,-10., -   2.,-   4., -   4.,    2.,     0.,    0., &
         17.,- 9., -  10.,   24.,     0.,    0.,     0.,    0., &
         1.,- 3., -   3.,-   1., -   2.,    5.,     0.,    0./
    DATA Z9/     1.,- 2., - 155.,-  52., -  78.,  193.,     7.,    0., &
         1.,- 1., -7208.,   59.,    56., 7067., -   1.,   17., &
         1.,  0., - 307.,-2582.,   227.,-  89.,    16.,    0., &
         1.,  1.,     8.,-  73.,    79.,    9.,     1.,   23., &
         2.,- 3.,    11.,   68.,   102.,-  17.,     0.,    0., &
         2.,- 2.,   136., 2728.,  4021.,- 203.,     0.,    0., &
         2.,- 1., - 537., 1518.,  1376.,  486.,    13.,  166., &
         2.,  0., -  22.,-  70., -   1.,-   8.,     0.,    0., &
         3.,- 4., -   5.,    2.,     3.,    8.,     0.,    0., &
         3.,- 3., - 162.,   27.,    43.,  278.,     0.,    0./
    DATA ZA/     3.,- 2.,    71.,  551.,   796.,- 104.,     6.,-   1., &
         3.,- 1., -  31.,  208.,   172.,   26.,     1.,   18., &
         4.,- 4., -   3.,-  16., -  29.,    5.,     0.,    0., &
         4.,- 3., -  43.,    9.,    13.,   73.,     0.,    0., &
         4.,- 2.,    17.,   78.,   110.,-  24.,     0.,    0., &
         4.,- 1., -   1.,   23.,    17.,    1.,     0.,    0., &
         5.,- 5.,     0.,    0., -   1.,-   3.,     0.,    0., &
         5.,- 4., -   1.,-   5., -  10.,    2.,     0.,    0., &
         5.,- 3., -   7.,    2.,     3.,   12.,     0.,    0., &
         5.,- 2.,     3.,    9.,    13.,-   4.,     0.,    0./
    DATA ZB/     1.,- 2., -   3.,   11.,    15.,    3.,     0.,    0., &
         1.,- 1., -  77.,  412.,   422.,   79.,     1.,    6., &
         1.,  0., -   3.,- 320.,     8.,-   1.,     0.,    0., &
         1.,  1.,     0.,-   8.,     8.,    0., -   1.,    6., &
         2.,- 3.,     0.,    0., -   3.,-   1.,     0.,    0., &
         2.,- 2.,    38.,- 101., - 152.,-  57.,     0.,    0., &
         2.,- 1.,    45.,- 103., - 103.,-  44.,    17., - 29., &
         2.,  0.,     2.,-  17.,     0.,    0.,     0.,    0., &
         3.,- 2.,     7.,-  20., -  30.,-  11.,     0.,    0., &
         3.,- 1.,     6.,-  16., -  16.,-   6.,     0.,    0./
    DATA ZC/     4.,- 2.,     1.,-   3., -   4.,-   1.,     0.,    0./

    CC(1,1) = +0.9999257180268403D0
    CC(1,2) = -0.0111782838886141D0
    CC(1,3) = -0.0048584357372842D0
    CC(2,1) = +0.0111782838782385D0
    CC(2,2) = +0.9999375206641666D0
    CC(2,3) = -0.0000271576311202D0
    CC(3,1) = +0.0048584357611567D0
    CC(3,2) = -0.0000271533600777D0
    CC(3,3) = +0.9999881973626738D0

    T  = (DJ - 2415020.D0)/36525.D0
    TP = (DJ - 2396758.D0)/365.25D0
    T18= (DJ - 2378496.D0)/36525.D0
    V = (DJ - 2451545.D0)/36525.D0
    V2 = V*V
    V3 = V*V2
    T50= (DJ - 2433282.423D0)/36525.D0
    EL( 2) =  0.16750401D-01 - T*(0.41879D-04 + 0.0503D-06*T)
    EL( 5) =  0.4908229466868880D+01 + (0.3000526416797348D-01 &
         +(0.7902463002085429D-05 +  0.5817764173314427D-07*T)*T)*T &
         -(5.066D3 + 5.880D3*T + 0.476D3*T*T)/STR
    EL( 9) =  0.4881627934111871D+01 + (0.6283319509909086D+03 &
         + 0.5279620987282842D-05*T)*T &
         +(0.33211D3 - 0.05299D3*T + 0.03959D3*T*T)/STR
    EL(10) =  0.6256583580497099D+01 + (0.6283019457267408D+03 &
         -(0.2617993877991491D-05 +  0.5817764173314427D-07*T)*T)*T &
         +(5.358D3 + 5.747D3*T + 0.516D3*T*T)/STR
    EL(14) =  0.5859485105530519D+01 + (0.8399709200339593D+04 &
         +(0.4619304753611657D-04 +  0.6544984694978728D-07*T18) &
         *T18)*T18
    EL(17) =  0.5807638839584459D+00 - (0.3375728238223387D+02 &
         -(0.3964806284113784D-04 +  0.3470781143063166D-07*T18) &
         *T18)*T18
    EL(21) =  0.3933938487925745D+01 + (0.7101833330913285D+02 &
         -(0.1752456013106637D-03 -  0.1773109075921903D-06*T18) &
         *T18)*T18
    EL(15) = EL(14) - EL(21)
    EL(16) = EL(14) - EL(17)
    D      = EL(14) - EL(9)
    ARG = D
    DL =    + 6466.D0*DSIN(ARG) + 13.D0*DSIN(3.*ARG)
    DR =    + 13384.D0*DCOS(ARG) +  30.D0*DCOS(3.*ARG)
    DBLARG = D + EL(15)
    ARG = DBLARG
    DL = DL +  177.D0*DSIN(ARG)
    DR = DR +  370.D0*DCOS(ARG)
    DBLARG = D - EL(15)
    ARG = DBLARG
    DL = DL -  425.D0*DSIN(ARG)
    DR = DR - 1332.D0*DCOS(ARG)
    DBLARG = 3.D0*D - EL(15)
    ARG = DBLARG
    DL = DL +   39.D0*DSIN(ARG)
    DR = DR +   80.D0*DCOS(ARG)
    DBLARG = D + EL(10)
    ARG = DBLARG
    DL = DL -   64.D0*DSIN(ARG)
    DR = DR -  140.D0*DCOS(ARG)
    DBLARG = D - EL(10)
    ARG = DBLARG
    DL = DL +  172.D0*DSIN(ARG)
    DR = DR + 360.D0*DCOS(ARG)
    ARG = D - EL(10) - EL(15)
    DL = DL - 13.D0*DSIN(ARG)
    DR = DR - 30.D0*DCOS(ARG)
    ARG = 2.D0*(EL(9)-EL(17))
    DL = DL - 13.D0*DSIN(ARG)
    DR = DR + 30.D0*DCOS(ARG)
    ARG = EL(16)
    DB =    + 577.D0*DSIN(ARG)
    DBLARG = EL(16) + EL(15)
    ARG = DBLARG
    DB = DB +  16.D0*DSIN(ARG)
    DBLARG = EL(16) - EL(15)
    ARG = DBLARG
    DB = DB -  47.D0*DSIN(ARG)
    DBLARG = EL(16) - 2.D0*(EL(9) - EL(17))
    ARG = DBLARG
    DB = DB +  21.D0*DSIN(ARG)
    ARG = ARG - EL(15)
    DB = DB + 5.D0*DSIN(ARG)
    ARG = EL(10) + EL(16)
    DB = DB + 5.D0*DSIN(ARG)
    ARG = EL(10) - EL(16)
    DB = DB + 5.D0*DSIN(ARG)
    GME = 0.43296383D+01 + .2608784648D+02*TP
    GV  = 0.19984020D+01 + .1021322923D+02*TP
    GMA = 0.19173489D+01 + .3340556174D+01*TP
    GJ  = 0.25836283D+01 + .5296346478D+00*TP
    GS  = 0.49692316D+01 + .2132432808D+00*TP
    GJ = GJ + 0.579904067D-02 * DSIN(5.D0*GS - 2.D0*GJ + 1.1719644977D0 - &
         0.397401726D-03 * TP)
    DG = 266.*DSIN(.555015D0 + 2.076942D0*T) &
         +6400.*DSIN(4.035027D0 + .3525565D0*T) &
         +(1882.D0-16.*T)*DSIN(.9990265D0 + 2.622706D0*T)
    EL(10) = DG/STR + EL(10)
    EL(12) =   DSIN(     EL(10)) * (6909794. - (17274. + 21.*T)*T) &
         + DSIN(2.D0*EL(10)) * (  72334. -    362.*T) &
         + DSIN(3.D0*EL(10)) * (   1050. -      8.*T) &
         + DSIN(4.D0*EL(10)) *       17.
    RO     =                          30594. -    152.*T &
         - DCOS(     EL(10)) * (7273841. - (18182. + 22.*T)*T) &
         - DCOS(2.D0*EL(10)) * (  91374. -    457.*T) &
         - DCOS(3.D0*EL(10)) * (   1446. -     11.*T) &
         - DCOS(4.D0*EL(10)) *       25.
    EL(11) = 10.D0**(RO*1.D-09)
    DO K=1,4
       DBLARG = Z(1,K)*GME + Z(2,K)*EL(10)
       ARG = DBLARG
       CS  = DCOS(ARG)
       SS  = DSIN(ARG)
       DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
       DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
    ENDDO
    DO K=5,44
       DBLARG = Z(1,K)*GV + Z(2,K)*EL(10)
       ARG = DBLARG
       CS  = DCOS(ARG)
       SS  = DSIN(ARG)
       DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
       DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
       DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
    ENDDO
    DO K=45,89
       DBLARG = Z(1,K)*GMA + Z(2,K)*EL(10)
       ARG = DBLARG
       CS  = DCOS(ARG)
       SS  = DSIN(ARG)
       DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
       DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
       DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
    ENDDO
    DO K=90,110
       DBLARG = Z(1,K)*GJ + Z(2,K)*EL(10)
       ARG = DBLARG
       CS  = DCOS(ARG)
       SS  = DSIN(ARG)
       DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
       DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
       DB  =(Z(7,K)*CS  + Z(8,K)*SS) +DB
    ENDDO
    DO K=111,121
       DBLARG = Z(1,K)*GS + Z(2,K)*EL(10)
       ARG = DBLARG
       CS  = DCOS(ARG)
       SS  = DSIN(ARG)
       DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
       DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
       DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
    ENDDO
    C(1) = EL(11)*10.D0**(DR*1.D-09)
    C(4) = (DL + DG + EL(12))/STR + EL(9)
    C(5) = DB/STR
    C(2) = C(4)*RTD
    C(3) = C(5)*RTD
    BS = (-0.03903D3 + 0.03450D3*T - 0.00077D3*T*T)/STR
    BC = (+0.03363D3 + 0.01058D3*T + 0.00286D3*T*T)/STR
    BR = -0.01914D3/STR
    C(5) = C(5) + BS*DSIN(C(4)) + BC*DCOS(C(4)) +BR*DSIN(C(4)-EL(17))
    RS = (+0.01015D3 + 0.00035D3*T + 0.00023D3*T*T)/STR
    RC = (-0.02652D3 - 0.00002D3*T + 0.00054D3*T*T)/STR
    C(1) = C(1) + RS*DSIN(EL(10)) - RC*DCOS(EL(10))
    T2 = T*T
    T3 = T2*T
    OBL = OBL0 - C1*T - C2*T2 + C3*T3
    OBL = OBL/RTD
    RIN(1) = C(1)*DCOS(C(4))*DCOS(C(5))
    RIN(2) = C(1)*(DCOS(C(5))*DSIN(C(4))*DCOS(OBL)-DSIN(C(5))*DSIN(OBL))
    RIN(3) = C(1)*(DCOS(C(5))*DSIN(C(4))*DSIN(OBL)+DSIN(C(5))*DCOS(OBL))
    RDPSC=4.848136811D-06
    DPTTY=365242.2D0
    SDL=DJ
    IF (DJ-TFIN .lt. 0.) goto 101
    IF (DJ-TFIN .eq. 0.) goto 101
    IF (DJ-TFIN .gt. 0.) goto 100
100 CONTINUE
    SDL=TFIN
101 CONTINUE
    DT=(SDL-2415020.5D0)/DPTTY
    DT2=DT**2
    DY=DABS(DJ-TFIN)/DPTTY
    DY2=DY**2
    DY3=DY**3
    AK=((23042.53D0+139.73D0*DT+.06D0*DT2)*DY+(30.23D0-.27D0*DT)*DY2+ &
         18.00D0*DY3)*RDPSC
    AW=((23042.53D0+139.73D0*DT+.06D0*DT2)*DY+(109.50D0+.39D0*DT)*DY2+ &
         18.32D0*DY3)*RDPSC
    BN=((20046.85D0-85.33D0*DT-.37D0*DT2)*DY+(-42.67D0-.37D0*DT)*DY2- &
         41.80D0*DY3)*RDPSC
    SINK=DSIN(AK)
    COSK=DCOS(AK)
    SINW=DSIN(AW)
    COSW=DCOS(AW)
    SINN=DSIN(BN)
    COSN=DCOS(BN)
    X(1) =-SINK*SINW   +COSK*COSW*COSN
    Y(1) =-COSK*SINW   -SINK*COSW*COSN
    W(1) =-COSW*SINN
    X(2) =SINK*COSW    +COSK*SINW*COSN
    Y(2) =COSK*COSW    -SINK*SINW*COSN
    W(2) =-SINW*SINN
    X(3) =COSK*SINN
    Y(3) =-SINK*SINN
    W(3) =COSN
    ROU(1) = 0.D0
    ROU(2) = 0.D0
    ROU(3) = 0.D0
    IF (DJ-TFIN .lt. 0.) goto 1
    IF (DJ-TFIN .eq. 0.) goto 1
    IF (DJ-TFIN .gt. 0.) goto 2
1   CONTINUE
    DO I=1,3
       ROU(I)=X(I)*RIN(1) + Y(I)*RIN(2)+W(I)*RIN(3)
    ENDDO
    GO TO 3
2   CONTINUE
    ROU(1)=      X(1)*RIN(1)   +X(2)*RIN(2)   +X(3)*RIN(3)
    ROU(2)=      Y(1)*RIN(1)   +Y(2)*RIN(2)   +Y(3)*RIN(3)
    ROU(3)=      W(1)*RIN(1)   +W(2)*RIN(2)   +W(3)*RIN(3)
3   CONTINUE
    RIH = DSQRT(ROU(1)*ROU(1)+ROU(2)*ROU(2)+ROU(3)*ROU(3))
    RIG = DATAN2(ROU(2),ROU(1))
    DIG = DATAN2(ROU(3),DSQRT(ROU(1)*ROU(1)+ROU(2)*ROU(2)))
    RIG = RIG + 0.0D3/STR
    ROU(1) = RIH*DCOS(RIG)*DCOS(DIG)
    ROU(2) = RIH*DSIN(RIG)*DCOS(DIG)
    DO I = 1,3
       ROUT(I) = CC(I,1)*ROU(1)+CC(I,2)*ROU(2)+CC(I,3)*ROU(3)
    ENDDO
    RETURN
  END subroutine newcomb

end module xvutils
