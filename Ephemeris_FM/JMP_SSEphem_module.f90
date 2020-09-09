module jmp_planetlib
implicit none
!**************************************************************************************************************************
! Fortran module preapared by F. Mignard for J.M. Petit in January 2019
! Object: Single autonomous module to compute relatively low accuracy planetary ephemeris around year 2000 +/- 200 years
! Built out of FM  math_astro.module and planet.module
! Main routine reference
! Reference:  Simon, J.L, Bretagnon, P., Chapront, J., &
!              Chapront-Touze, M., Francou, G., and Laighskar, J., &
!              Astron. Astrophys. 282, 663 (1994).
! Version     10 January  2019
!  
! Only one call is useful for the user when the module is linked to the main (and with 'use jmp_planetlib' in the calling routines)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! call ephem_planet_simon(datjd, ipla, iframe, rr, vv)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Ephemeris of the planets and sun in the BCRS frame (Equatorial or eclipic) with the origin at the SS barycenter.
!     Source : Approximate osculating elements of Simon et al.
!     accuracy    5 to 10 "  ,  5m/s
!     datjd     : date in JD. (24xxxxx.yyyy)
!     iplan     : index of the planet  0 : Sun, 1 : Mercury, ... 8 : Neptune, 9 : EMB
!     iframe    :  0   =  Equatorial frame, 1 = ecliptic frame (ICRS, ie close to J2000)
!     flvit     : .true. ==> velocity vector computed, .false. ==> vv = 0
!     rr        : position vector   in km.  rr(3)
!     vv        : velocity vector   in km/s vv(3)
!
! The routine
! call  planet_helio_simon(iref, iplan, datjd, xlong, xlat, dist, pos, vit)
! gives the ephemeris in heliocentric frame with respect to J2000 ecliptic (iref =0) or mean ecliptic & equinox of date (iref=1)
!**************************************************************************************************************************
!
!*********************************************************
! Constantes

         integer,       parameter :: dp        = selected_real_kind(15,60)      ! real(8) in 64 bits
         real(kind=dp), parameter :: xpi       = 3.141592653589793238462643d0
         real(kind=dp), parameter :: deupi     = 6.283185307179586476925286d0
         real(kind=dp), parameter :: round     = 360d0
         real(kind=dp), parameter :: degrad    = xpi/180d0                      !  1.74532925199433d-2   :: degree to radian
         real(kind=dp), parameter :: raddeg    = 180d0/xpi                      !  57.2957795130823d0    :: radian to degree
         real(kind=dp), parameter :: xjd2000   = 2451545d0                      ! 1.0 Jan 2000 12h = J2000.0
         real(kind=dp), parameter :: xjd2010   = 2455197.5d0                    ! 1.0 Jan 2010 0h  = J2010.0.
         real(kind=dp), parameter :: days      = 86400d0                        ! day in seconds
         real(kind=dp), parameter :: yearj     = 365.25d0                       ! julian year    in days
         real(kind=dp), parameter :: century   = 36525d0                        ! julian century in days
         real(kind=dp), parameter :: asdeg     = 1d0/3600d0                     !  2.77777777777777d-4   :: second to degrees
         real(kind=dp), parameter :: xau       = 149597870.700d0                !  AU in km (TDB) 
         real(kind=dp), parameter :: vlight    = 299792.458d0                   ! velocity of light in km s^-1)
         real(kind=dp), parameter :: audtokms  = xau/days                       !  1731.45683670139d0    :: au/day to km/s
         character(10), dimension(0:9), parameter :: ssname =&                             ! Names of the planets
                                (/' Sun      ',' Mercury  ',' Venus    ',' Earth    ',' Mars     ', ' Jupiter  ',&
                                 ' Saturn   ', ' Uranus   ', ' Neptune  ', ' EMB      '/)

contains
!********************************************************************
   subroutine ephem_planet_simon(datjd, ipla, iframe, rr, vv)
 !********************************************************************
      !     Ephemeris of the planets and sun in the ICRS frame with the origin at the SS barycenter.
      !     Source : Approximate osculating elements of Simon et al.
      !     accuracy    5 to 10 "  ,  5m/s
      !     datjd     : date in JD
      !     iplan     : index of the planet  0 : Sun, 1 : Mercury, ... 8 : Neptune, 9 : EMB
      !     iframe    :  0   =  Equatorial frame, 1 = ecliptic frame
      !     flvit     : .true. ==> velocity vector computed, .false. ==> vv = 0
      !     rr        : position vector   in km
      !     vv        : velocity vector   in km/s
	  !	  !
	  !      Author  : F. Mignard  September 2015 from core routine provided by J.L. Simon
	  !
	  !      Source  : 
	  !      The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront, &
	  !      M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
	  !      Longitudes, Paris, France).  From comparisons with JPL
	  !      ephemeris DE102, they quote the following maximum errors
	  !      over the interval 1800-2050:
 
      !
      !*******************************************************************
      !     Comparisons against DE200 over the interval 1800-2100 gave the
      !     following maximum absolute differences.  (The results using
      !     DE406 were essentially the same.)
      !
      !                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)
      !
      !        Mercury        7            1            500       0.7
      !        Venus          7            1           1100       0.9
      !        EMB            9            1           1300       1.0
      !        Mars          26            1           9000       2.5
      !        Jupiter       78            6          82000       8.2
      !        Saturn        87           14         263000      24.6
      !        Uranus        86            7         661000      27.4
      !        Neptune       11            2         248000      21.4
      !
      !
      !********************************************************************
      implicit none

      real(kind = dp), parameter :: xmse = 1.0d0/82.30056d0 ! mass moon/(earth + moon)

      integer, parameter :: ieqecl = 1 ! Precession in ecliptic frame
      integer, parameter :: iopt = 1 ! Geometric position for the Moon

      integer, parameter :: iprec = 0 ! low precision for the Moon

      real(kind = dp) rr(3), vv(3), xv(6)
      real(kind = dp) posm(3), vitm(3)
      real(kind = dp) xp(3), xpp(3)
      real(kind = dp) xyz(3), xyzp(3)
      real(kind = dp) xlong, xlat, dist, ra, dec, rho
      real(kind = dp) datjd
      integer iplan, iframe, ierr, ipla

      iplan = ipla
      if (ipla == 9) iplan = 3 ! for EMB
      !
      !     position vector
      !
      if (iplan < 4) then
         call baryc_m(datjd, xyz, xyzp) !sun wrt. SS barycenter in au Ecliptic J2000  medium accuracy
      else
         call baryc_m(datjd, xyz, xyzp) !sun wrt. SS barycenter in au Ecliptic J2000 low accuracy
      endif

      rr = xyz * xau !sun wrt. SS barycenter in km
      vv = xyzp * xau/days !sun wrt. SS barycenter in km/s

      if (iplan > 0) then ! all the planets
         call planetap_sofa(datjd, iplan, xlong, xlat, dist, xv, ierr) ! au, au/day heliocentric of date
         xp = xv(1:3) * xau ! km
         xpp = xv(4:6) * xau/days ! km/s

         rr = rr + xp ! pos BCRS J2000 in km
         vv = vv + xpp ! vit in km/s
      endif

      if (ipla == 3) then ! EMB to Earth

         call Moon(datjd, iprec, iopt, xlong, xlat, rho, ra, dec, posm) ! Moon mean equinox and ecliptic, km
         call sphecar(rho, xlong, xlat, posm) ! posm in the procedure is in equatorial frame, now in ecliptic
         call preces(ieqecl, datjd, xjd2000, posm, posm) ! J2000 frame

         rr = rr - posm * xmse ! pos in km BCRS J2000

         call vit_moon(datjd, iframe, vitm)
         call preces(ieqecl, datjd, xjd2000, vitm, vitm) ! J2000 frame

         vv = vv - vitm * xmse ! km/s omegaXr is neglected for the EMB( < 1mm/s)

      endif

      if (iframe == 0) then
         call ecltoequ(xjd2000, rr, rr) ! passage to equatorial ICRS frame
         call ecltoequ(xjd2000, vv, vv) ! passage to equatorial ICRS frame
      endif

      return
   end subroutine ephem_planet_simon
!********************************************************************
subroutine planetap_sofa(datjd, ipla, xlong, xlat, dist, xv, ierr)
   !********************************************************************
   !
   !     interface to use IAU_PLAN94 of SOFA with the same call list as Planetap
   !     in this library.
   !
   !     Approximate heliocentric position and velocity of a nominated major
   !     planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
   !     Neptune (but not Pluto nor the Earth itself).
   !     (heliocentric, J2000, AU, AU/d) Ecliptic frame (FM 08/04)!
   !
   !     Input
   !           datjd    Date in julian days
   !           ipla     Index of the planet 1 : Mercury , 3 : EMB, .. 8 : Neptune
   !           xlong    Ecliptic longitude in ecliptic J2000 in degrees - heliocentric
   !           xlat     Ecliptic latitude  in ecliptic J2000 in degrees - heliocentric
   !           dist     Radius vector in au - heliocentric
   !           xv       array xv(6) with xv(1:3) for position and xv(4:6) for velocity AU, AU/day
   !           ierr     status: -1 = illegal ipla (outside 1-8)
   !                             0 = OK
   !                            +1 = warning: date outside 1000-3000 AD
   !                            +2 = warning: solution failed to converge
   !
   !     From Sofa coding :
   !     Comparisons against DE200 over the interval 1800-2100 gave the
   !     following maximum absolute differences.  (The results using
   !     DE406 were essentially the same.)
   !
   !                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)
   !
   !        Mercury        7            1            500       0.7
   !        Venus          7            1           1100       0.9
   !        EMB            9            1           1300       1.0
   !        Mars          26            1           9000       2.5
   !        Jupiter       78            6          82000       8.2
   !        Saturn        87           14         263000      24.6
   !        Uranus        86            7         661000      27.4
   !        Neptune       11            2         248000      21.4
   !
   !
   !********************************************************************
   implicit none

   real(kind = dp) datjd, xlong, xlat, dist
   real(kind = dp) xv(6)

   real(kind = dp) epoch2
   real(kind = dp) pv(3, 2)

   integer ipla, ierr

   epoch2 = 0d0

   call plan_simon(datjd, epoch2, ipla, pv, ierr)

   call carsphe(pv(:, 1), dist, xlong, xlat)

   xv(1:3) = pv(:, 1)
   xv(4:6) = pv(:, 2)

   return
end subroutine planetap_sofa

!********************************************************************
subroutine plan_simon(epoch1, epoch2, np, pv, j)
!********************************************************************
   !
   !  Approximate heliocentric position and velocity of a nominated major
   !  planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
   !  Neptune (but not Pluto nor the Earth itself).
   !  SOFA implementation (replaces FM version of 2002)
   !
   !  Given:
   !     EPOCH1   d       TDB epoch part A (Note 1)
   !     EPOCH2   d       TDB epoch part B (Note 1)
   !     NP       i       planet (1=Mercury, 2=Venus, 3=EMB ... 8=Neptune)
   !
   !  Returned:
   !     PV       d(3,2)  planet pos,vel (heliocentric, J2000, AU, AU/d) Ecliptic frame (FM 08/04)
   !     J        i       status: -1 = illegal NP (outside 1-8)
   !                               0 = OK
   !                              +1 = warning: date outside 1000-3000 AD
   !                              +2 = warning: solution failed to converge

   !  Notes
   !
   !  1) The epoch EPOCH1+EPOCH2 is in the TDB timescale and is a Julian
   !     Date, apportioned in any convenient way between the two arguments.
   !     For example, JD(TDB)=2450123.7 could be expressed in any of these
   !     ways, among others:
   !
   !            EPOCH1        EPOCH2
   !
   !         2450123.7D0        0D0        (JD method)
   !         2451545D0       -1421.3D0     (J2000 method)
   !         2400000.5D0     50123.2D0     (MJD method)
   !         2450123.5D0       0.2D0       (date & time method)
   !
   !     The JD method is the most natural and convenient to use in
   !     cases where the loss of several decimal digits of resolution
   !     is acceptable.  The J2000 method is best matched to the way
   !     the argument is handled internally and will deliver the
   !     optimum resolution.  The MJD method and the date & time methods
   !     are both good compromises between resolution and convenience.
   !     The limited accuracy of the present algorithm is such that any
   !     of the methods is satisfactory.
   !
   !  2) If an NP value outside the range 1-8 is supplied, an error
   !     status (J = -1) is returned and the PV vector set to zeroes.
   !
   !  3) For NP=3 the result is for the Earth-Moon Barycenter.  To
   !     obtain the heliocentric position and velocity of the Earth, &
   !     use instead the SOFA routine iau_EPV00.
   !
   !  4) On successful return, the array PV contains the following:
   !
   !        PV(1,1)  x       }
   !        PV(2,1)  y       } heliocentric position, AU
   !        PV(3,1)  z       }
   !
   !        PV(1,2)  xdot    }
   !        PV(2,2)  ydot    } heliocentric velocity, AU/d
   !        PV(3,2)  zdot    }
   !
   !     The reference frame is equatorial and is with respect to the
   !     mean equator and equinox of epoch J2000.
   !
   !  5) The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront, &
   !     M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
   !     Longitudes, Paris, France).  From comparisons with JPL
   !     ephemeris DE102, they quote the following maximum errors
   !     over the interval 1800-2050:
   !
   !                     L (arcsec)    B (arcsec)      R (km)
   !
   !        Mercury          4             1             300
   !        Venus            5             1             800
   !        EMB              6             1            1000
   !        Mars            17             1            7700
   !        Jupiter         71             5           76000
   !        Saturn          81            13          267000
   !        Uranus          86             7          712000
   !        Neptune         11             1          253000
   !
   !     Over the interval 1000-3000, they report that the accuracy is no
   !     worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
   !     accuracy declines.
   !
   !     Comparisons of the present routine with the JPL DE200 ephemeris
   !     give the following RMS errors over the interval 1960-2025:
   !
   !                      position (km)     velocity (m/s)
   !
   !        Mercury            334               0.437
   !        Venus             1060               0.855
   !        EMB               2010               0.815
   !        Mars              7690               1.98
   !        Jupiter          71700               7.70
   !        Saturn          199000              19.4
   !        Uranus          564000              16.4
   !        Neptune         158000              14.4
   !
   !     Comparisons against DE200 over the interval 1800-2100 gave the
   !     following maximum absolute differences.  (The results using
   !     DE406 were essentially the same.)
   !
   !                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)
   !
   !        Mercury        7            1            500       0.7
   !        Venus          7            1           1100       0.9
   !        EMB            9            1           1300       1.0
   !        Mars          26            1           9000       2.5
   !        Jupiter       78            6          82000       8.2
   !        Saturn        87           14         263000      24.6
   !        Uranus        86            7         661000      27.4
   !        Neptune       11            2         248000      21.4
   !
   !  6) The present SOFA re-implementation of the original Simon et al.
   !     Fortran code differs from the original in the following respects:
   !
   !       *  The date is supplied in two parts.
   !
   !
   !       *  More is done in-line:  there are fewer calls to other
   !          routines.
   !
   !       *  Different error/warning status values are used.
   !
   !       *  A different Kepler's-equation-solver is used (avoiding
   !          use of COMPLEX*16).
   !
   !       *  Polynomials in T are nested to minimize rounding errors.
   !
   !       *  Explicit double-precision constants are used to avoid
   !          mixed-mode expressions.
   !
   !       *  There are other, cosmetic, changes to comply with SOFA
   !          style conventions.
   !
   !     None of the above changes affects the result significantly.
   !
   !  7) The returned status, J, indicates the most serious condition
   !     encountered during execution of the routine.  Illegal NP is
   !     considered the most serious, overriding failure to converge, &
   !     which in turn takes precedence over the remote epoch warning.
   !
   !  Called:
   !     iau_ANP     normalize radians to range -pi to +pi
   !
   !  Reference:  Simon, J.L, Bretagnon, P., Chapront, J., &
   !              Chapront-Touze, M., Francou, G., and Laskar, J., &
   !              Astron. Astrophys. 282, 663 (1994).
   !
   !  This revision:  2001 May 24
   !
   !  Copyright (C) 2003 IAU SOFA Review Board.  See notes at end.
   !
   !-----------------------------------------------------------------------
   !********************************************************************
   implicit none

   real(kind=dp), intent(in)      :: epoch1, epoch2
   integer,       intent(in)      :: np
   real(kind=dp)                  :: pv(3, 2)
   integer                        :: j

   !  maximum number of iterations allowed to solve kepler's equation
   integer , parameter            :: kmax = 10

   !  2pi
   real(kind=dp), parameter       :: d2pi = 6.283185307179586476925287d0

   !  arc seconds to radians
   real(kind=dp), parameter       :: das2r = 4.848136811095359935899141d-6

   !  reference epoch (j2000), jd
   real(kind=dp), parameter       :: dj0 = 2451545d0

   !  days per julian millennium
   real(kind=dp), parameter       :: djm = 365250d0

   !  sin and cos of j2000 mean obliquity (iau 1976)
   real(kind=dp), parameter       ::  sineps = 0.3977771559319137d0, coseps = 0.9174820620691818d0

   !  gaussian constant
   real(kind=dp), parameter       :: gk = 0.017202098950d0

   integer                        :: jstat, k
   real(kind=dp)                  :: amas(8), a(3, 8), dlm(3, 8), e(3, 8)
   real(kind=dp)                  :: pi(3, 8), dinc(3, 8), omega(3, 8)
   real(kind=dp)                  :: kp(9, 8), ca(9, 8), sa(9, 8)
   real(kind=dp)                  :: kq(10, 8), cl(10, 8), sl(10, 8)
   real(kind=dp)                  :: t, da, dl, de, dpp, di, dom, dmu, arga, argl, am
   real(kind=dp)                  :: ae, dae, ae2, at, r, v, si2, xq, xp, tl, xsw
   real(kind=dp)                  :: xcw, xm2, xf, ci2, xms, xmc, xpxq2, x, y, z

   !      double precision anpm

   !  planetary inverse masses
   data amas / 6023600d0, 408523.5d0, 328900.5d0, 3098710d0, 1047.355d0, 3498.5d0, 22869d0, 19314d0 /

   !
   !  Tables giving the mean Keplerian elements, limited to T**2 terms:
   !
   !         A       semi-major axis (AU)
   !         DLM     mean longitude (degree and arcsecond)
   !         E       eccentricity
   !         PI      longitude of the perihelion (degree and arcsecond)
   !         DINC    inclination (degree and arcsecond)
   !         OMEGA   longitude of the ascending node (degree and arcsecond)
   !
   DATA A /&
   0.3870983098D0, 0D0, 0D0, &
   0.7233298200D0, 0D0, 0D0, &
   1.0000010178D0, 0D0, 0D0, &
   1.5236793419D0, 3d-10, 0D0, &
   5.2026032092D0, 19132d-10, -39d-10, &
   9.5549091915D0, -0.0000213896D0, 444d-10, &
   19.2184460618D0, -3716d-10, 979d-10, &
   30.1103868694D0, -16635d-10, 686d-10 /
   !
   DATA DLM /&
   252.25090552D0, 5381016286.88982D0, -1.92789D0, &
   181.97980085D0, 2106641364.33548D0, 0.59381D0, &
   100.46645683D0, 1295977422.83429D0, -2.04411D0, &
   355.43299958D0, 689050774.93988D0, 0.94264D0, &
   34.35151874D0, 109256603.77991D0, -30.60378D0, &
   50.07744430D0, 43996098.55732D0, 75.61614D0, &
   314.05500511D0, 15424811.93933D0, -1.75083D0, &
   304.34866548D0, 7865503.20744D0, 0.21103D0 /
   !
   DATA E /&
   0.2056317526D0, 0.0002040653D0, -28349d-10, &
   0.0067719164D0, -0.0004776521D0, 98127d-10, &
   0.0167086342D0, -0.0004203654D0, -0.0000126734D0, &
   0.0934006477D0, 0.0009048438D0, -80641d-10, &
   0.0484979255D0, 0.0016322542D0, -0.0000471366D0, &
   0.0555481426D0, -0.0034664062D0, -0.0000643639D0, &
   0.0463812221D0, -0.0002729293D0, 0.0000078913D0, &
   0.0094557470D0, 0.0000603263D0, 0D0 /
   !
   DATA PI /&
   77.45611904D0, 5719.11590D0, -4.83016D0, &
   131.56370300D0, 175.48640D0, -498.48184D0, &
   102.93734808D0, 11612.35290D0, 53.27577D0, &
   336.06023395D0, 15980.45908D0, -62.32800D0, &
   14.33120687D0, 7758.75163D0, 259.95938D0, &
   93.05723748D0, 20395.49439D0, 190.25952D0, &
   173.00529106D0, 3215.56238D0, -34.09288D0, &
   48.12027554D0, 1050.71912D0, 27.39717D0 /
   !
   DATA DINC /&
   7.00498625D0, -214.25629D0, 0.28977D0, &
   3.39466189D0, -30.84437D0, -11.67836D0, &
   0D0, 469.97289D0, -3.35053D0, &
   1.84972648D0, -293.31722D0, -8.11830D0, &
   1.30326698D0, -71.55890D0, 11.95297D0, &
   2.48887878D0, 91.85195D0, -17.66225D0, &
   0.77319689D0, -60.72723D0, 1.25759D0, &
   1.76995259D0, 8.12333D0, 0.08135D0 /
   !
   DATA OMEGA /&
   48.33089304D0, -4515.21727D0, -31.79892D0, &
   76.67992019D0, -10008.48154D0, -51.32614D0, &
   174.87317577D0, -8679.27034D0, 15.34191D0, &
   49.55809321D0, -10620.90088D0, -230.57416D0, &
   100.46440702D0, 6362.03561D0, 326.52178D0, &
   113.66550252D0, -9240.19942D0, -66.23743D0, &
   74.00595701D0, 2669.15033D0, 145.93964D0, &
   131.78405702D0, -221.94322D0, -0.78728D0 /
   !
   !  Tables for trigonometric terms to be added to the mean elements
   !  of the semi-major axes.
   !
   DATA KP /&
   69613, 75645, 88306, 59899, 15746, 71087, 142173, 3086, 0, &
   21863, 32794, 26934, 10931, 26250, 43725, 53867, 28939, 0, &
   16002, 21863, 32004, 10931, 14529, 16368, 15318, 32794, 0, &
   6345, 7818, 15636, 7077, 8184, 14163, 1107, 4872, 0, &
   1760, 1454, 1167, 880, 287, 2640, 19, 2047, 1454, &
   574, 0, 880, 287, 19, 1760, 1167, 306, 574, &
   204, 0, 177, 1265, 4, 385, 200, 208, 204, &
   0, 102, 106, 4, 98, 1367, 487, 204, 0 /
   !
   DATA CA /&
   4, -13, 11, -9, -9, -3, -1, 4, 0, &
   -156, 59, -42, 6, 19, -20, -10, -12, 0, &
   64, -152, 62, -8, 32, -41, 19, -11, 0, &
   124, 621, -145, 208, 54, -57, 30, 15, 0, &
   -23437, -2634, 6601, 6259, -1507, -1821, 2620, -2115, -1489, &
   62911, -119919, 79336, 17814, -24241, 12068, 8306, -4893, 8902, &
   389061, -262125, -44088, 8387, -22976, -2093, -615, -9720, 6633, &
   -412235, -157046, -31430, 37817, -9740, -13, -7449, 9644, 0 /
   !
   DATA SA /&
   -29, -1, 9, 6, -6, 5, 4, 0, 0, &
   -48, -125, -26, -37, 18, -13, -20, -2, 0, &
   -150, -46, 68, 54, 14, 24, -28, 22, 0, &
   -621, 532, -694, -20, 192, -94, 71, -73, 0, &
   -14614, -19828, -5869, 1881, -4372, -2255, 782, 930, 913, &
   139737, 0, 24667, 51123, -5102, 7429, -4095, -1976, -9566, &
   -138081, 0, 37205, -49039, -41901, -33872, -27037, -12474, 18797, &
   0, 28492, 133236, 69654, 52322, -49577, -26430, -3593, 0 /
   !
   !  Tables giving the trigonometric terms to be added to the mean
   !  elements of the mean longitudes.
   !
   DATA KQ /&
   3086, 15746, 69613, 59899, 75645, 88306, 12661, 2658, 0, 0, &
   21863, 32794, 10931, 73, 4387, 26934, 1473, 2157, 0, 0, &
   10, 16002, 21863, 10931, 1473, 32004, 4387, 73, 0, 0, &
   10, 6345, 7818, 1107, 15636, 7077, 8184, 532, 10, 0, &
   19, 1760, 1454, 287, 1167, 880, 574, 2640, 19, 1454, &
   19, 574, 287, 306, 1760, 12, 31, 38, 19, 574, &
   4, 204, 177, 8, 31, 200, 1265, 102, 4, 204, &
   4, 102, 106, 8, 98, 1367, 487, 204, 4, 102 /
   !
   DATA CL /&
   21, -95, -157, 41, -5, 42, 23, 30, 0, 0, &
   -160, -313, -235, 60, -74, -76, -27, 34, 0, 0, &
   -325, -322, -79, 232, -52, 97, 55, -41, 0, 0, &
   2268, -979, 802, 602, -668, -33, 345, 201, -55, 0, &
   7610, -4997, -7689, -5841, -2617, 1115, -748, -607, 6074, 354, &
   -18549, 30125, 20012, -730, 824, 23, 1289, -352, -14767, -2062, &
   -135245, -14594, 4197, -4030, -5630, -2898, 2540, -306, 2939, 1986, &
   89948, 2103, 8963, 2695, 3682, 1648, 866, -154, -1963, -283 /
   !
   DATA SL /&
   -342, 136, -23, 62, 66, -52, -33, 17, 0, 0, &
   524, -149, -35, 117, 151, 122, -71, -62, 0, 0, &
   -105, -137, 258, 35, -116, -88, -112, -80, 0, 0, &
   854, -205, -936, -240, 140, -341, -97, -232, 536, 0, &
   -56980, 8016, 1012, 1448, -3024, -3710, 318, 503, 3767, 577, &
   138606, -13478, -4964, 1441, -1319, -1482, 427, 1236, -9167, -1918, &
   71234, -41116, 5334, -4935, -1848, 66, 434, -1748, 3780, -701, &
   -47645, 11647, 2166, 3194, 679, 0, -244, -419, -2531, 48 /

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


   !  Validate the planet number.
   if (np .lt. 1 .or. np .gt. 8) then
      jstat = -1
      pv = 0d0          !     reset the result in case of failure.
   else            
      t = ((epoch1 - dj0) + epoch2) / djm     !     time: julian millennia since j2000.
      if (abs(t) .le. 1d0) then
         jstat = 0
      else
         jstat = 1
      end if
!
!     compute the mean elements.
! 
      da = a(1, np) +(a(2, np) + a(3, np) * t) * t
      dl = (3600d0 * dlm(1, np) + (dlm(2, np) +  dlm(3, np) * t) * t) * das2r
      de = e(1, np) + (e(2, np) +   e(3, np) * t) * t
      dpp = xmodpi((3600d0 * pi(1, np) +  (pi(2, np) + pi(3, np) * t) * t) * das2r, d2pi)
      di = (3600d0 * dinc(1, np) +   (dinc(2, np) + dinc(3, np) * t) * t) * das2r
      dom = xmodpi((3600d0 * omega(1, np) +(omega(2, np) +omega(3, np) * t) * t) * das2r,d2pi)
!
!     apply the trigonometric terms
!
      dmu = 0.35953620d0 * t !  0.35953620d0  = (n5-n6)/880 in rad/millenium
      do  k = 1, 8
         arga = kp(k, np) * dmu
         argl = kq(k, np) * dmu
         da = da + (ca(k, np) * cos(arga) +  sa(k, np) * sin(arga)) * 1d-7
         dl = dl + (cl(k, np) * cos(argl) +   sl(k, np) * sin(argl)) * 1d-7
      enddo
      
      arga = kp(9, np) * dmu
      da = da + t * (ca(9, np) * cos(arga) +  sa(9, np) * sin(arga)) * 1d-7
      do  k = 9, 10
         argl = kq(k, np) * dmu
         dl = dl + t * (cl(k, np) * cos(argl) + sl(k, np) * sin(argl)) * 1d-7
      enddo
      dl = modulo(dl, d2pi)  
!      
!     iterative solution of kepler's equation to get eccentric anomaly.
!
      am = dl - dpp
      ae = am + de * sin(am)
      k = 0
            
      do
         dae = (am - ae + de * sin(ae)) / (1d0 - de * cos(ae))
         ae = ae + dae
         k = k + 1
         if (k .ge. kmax) jstat = 2
         if ((k .gt. kmax) .or. (abs(dae) .gt. 1d-12)) exit
      enddo  
      
  !  true anomaly.
      ae2 = ae / 2d0
      at = 2d0 * atan2(sqrt((1d0 + de)/(1d0 - de)) * sin(ae2),  cos(ae2))

  !  distance (au) and speed (radians per day).
      r = da * (1d0 - de * cos(ae))
      v = gk * sqrt((1d0 + 1d0/amas(np)) / (da * da * da))

      si2 = sin(di/2d0)
      xq = si2 * cos(dom)
      xp = si2 * sin(dom)
      tl = at + dpp
      xsw = sin(tl)
      xcw = cos(tl)
      xm2 = 2d0 * (xp * xcw - xq * xsw)
      xf = da / sqrt(1d0 - de * de)
      ci2 = cos(di/2d0)
      xms = (de * sin(dpp) + xsw) * xf
      xmc = (de * cos(dpp) + xcw) * xf
      xpxq2 = 2d0 * xp * xq

      !     position (j2000 ecliptic x,y,z in au).
      x = r * (xcw - xm2 * xp)
      y = r * (xsw + xm2 * xq)
      z = r * (-xm2 * ci2)

 !  rotate to equatorial.  deleted fm 08/04 for backward compatibility with previous version
 !      pv(1,1) = x
 !      pv(2,1) = y*coseps - z*sineps
 !      pv(3,1) = y*sineps + z*coseps

 !  position (j2000 ecliptic x,y,z in au).
      pv(1, 1) = x
      pv(2, 1) = y
      pv(3, 1) = z

 !      pos(1) = x
 !      pos(2) = y
 !      pos(3) = z

      !     velocity (j2000 ecliptic xdot,ydot,zdot in au/d).
      x = v * ((-1d0 + 2d0 * xp * xp) * xms + xpxq2 * xmc)
      y = v * ((1d0 - 2d0 * xq * xq) * xmc - xpxq2 * xms)
      z = v * (2d0 * ci2 * (xp * xms + xq * xmc))

 !  rotate to equatorial. deleted fm 08/04 for backward compatibility with previous version
 !      pv(1,2) = x
 !      pv(2,2) = y*coseps - z*sineps
 !      pv(3,2) = y*sineps + z*coseps

 !  velocity (j2000 ecliptic xdot,ydot,zdot in au/d).
      pv(1, 2) = x
      pv(2, 2) = y
      pv(3, 2) = z

!       vit(1) = x
!       vit(2) = y
!       vit(3) = z
         
   endif

!  return the status.
   j = jstat

   end subroutine plan_simon

 

!*****************************************************************
!*****************************************************************
subroutine Moon(xjd, iprec, iopt, xlong, xlat, rho, ra, dec, pos)
   !*****************************************************************
   !  Ephemeris of the Moon
   !  Based on ELP-2000 updated with LLR data
   !  Truncated version for historical purposes
   !  Source : ELP2000-85, a semi-analytical lunar ephemeris for
   !  for historical times, Chapront Touze M., Chapront J., 1988,  A&A, 190, 342.
   !
   !  Updated for the mean arguments from Chapront et al. A&A  2002, 387, 700-709 by FM in December 2005
   !  Precession constants updated by FM in December 2015 with updated IAU values
   !
   !  ICRS x-axis instead of equinox F. Mignard December 2005
   !
   !  Higher accuracy (1 mas) version with ELP-2000-85 in module lunar_ephem and routine 
   !  Mond(xjd, accur, iframe, icrs, xlong, xlat, rho, pos, vel)
   !
   !
   !  F. Mignard  December 2005
   !
   !  INPUT
   !   xjd      : date in julian days
   !   iprec    : accuracy level     0 : 5",10km , 1: 1.0", 2km,  2: 0.05", 0.1km (more terms computed)
   !   iopt     : 1 : geometric, mean equinox, ecliptic, equator of date
   !            : 2 : apparent,  true equinox and equator of date, true equinox, mean ecliptic of date.
   !
   !  OUTPUT
   !   xlong    : ecliptic longitude          deg
   !   xlat     : ecliptic latitude           deg
   !   rho      : geocentric distance         km
   !   ra       : right ascension             deg
   !   dec      : declination                 deg
   !   pos      : equatorial position vector  km
   !
   !
   !   Precision level determined by F. Mignard and set by truncating each series
   !   at the applicable level ( amplitude ~ accuracy/10 and t = 10 centuries for mixed terms)
   !   Coorresponding series maximum indexes are given in parameter.
   !
   !************************************************************

   implicit none

   real(kind = dp), intent(in) :: xjd
   integer, intent(in) :: iprec
   integer, intent(in) :: iopt
   real(kind = dp), intent(out) :: xlong, xlat, rho, ra, dec
   real(kind = dp), intent(out) :: pos(3)

   !
   ! Arrays with the number of terms in every series
   !

   integer, parameter :: nmaxcent = 218
   integer, dimension(3), parameter :: nmcent = (/218, 188, 155/)
   integer, dimension(3), parameter :: nmpla0 = (/244, 64, 115/)
   integer, dimension(3), parameter :: nmpla1 = (/154, 64, 69/)
   integer, dimension(3), parameter :: nmpla2 = (/ 25, 12, 19/)

   !
   ! Arrays giving the number of terms selected to reach the accuracy determined by iprec (computed for t = 10 centuries)
   !

   integer, dimension(3, 0:2), parameter :: idmain = reshape([70, 46, 42, 130, 95, 89, 218, 188, 155], [3, 3])
   integer, dimension(3, 0:2), parameter :: idpla0 = reshape([4, 3, 0, 50, 10, 11, 244, 64, 115], [3, 3])
   integer, dimension(3, 0:2), parameter :: idpla1 = reshape([8, 2, 5, 25, 15, 15, 154, 64, 69], [3, 3])
   integer, dimension(3, 0:2), parameter :: idpla2 = reshape([0, 0, 0, 2, 2, 2, 25, 12, 19], [3, 3])

   !
   ! Delaunay arguments Chapront et al. A&A  2002, 387, 700-709 D                l'                  l               F
   !
   real(kind = dp), dimension(4), parameter :: delau0 = (/ 297.85019167d0, 357.52910389d0, 134.96339622d0, 93.272097694d0/)
   real(kind = dp), dimension(4), parameter :: delau1 = (/ 1602961601.0312d0, 129596581.0733d0, 1717915923.0024d0, 1739527263.2179d0/)
   real(kind = dp), dimension(4), parameter :: delau2 = (/-6.8498d0, -0.5529d0, 31.3939d0, -13.2293d0/)
   real(kind = dp), dimension(4), parameter :: delau3 = (/ 0.006595d0, 0.000147d0, 0.051651d0, -0.001021d0/)
   real(kind = dp), dimension(4), parameter :: delau4 = (/-0.00003184d0, 0.00000015d0, -0.00024470d0, 0.00000417d0/)

   !
   ! Mean longitude of the moon equinox J2000 (gamma ICRS)   t^0            t                t^2         t^3          t^4
   ! Chapront et al. A&A  2002, 387, 700-709

   real(kind = dp), dimension(0:4), parameter :: xlmoon = (/785939.8782d0, 1732559343.332d0, -6.870d0, +0.006604d0, -0.00003169d0/)

   !
   ! Precession in ecliptic in " and cy                      t             t^2        t^3        t^4
   ! Second line : modified by FM 21/09/2015 with IAU 2006 values

!   real(kind = dp), dimension(1:7), parameter :: cpreces = (/5028.7275d0, 1.105394d0, 0.000076d0, -0.2353d-4, -0.00181d-5, 0.000175d-6, 0.000013d-7/)
   real(kind = dp), dimension(1:7), parameter :: cpreces = (/5028.796195d0, 110.54348d-2, 0.07964d-3, -0.23857d-4, -0.00383d-5, 0.000175d-6, 0.000013d-7/)


   !
   ! Planetary mean longitude
   !
   real(kind = dp), dimension(8), parameter :: planet0 = (/252.250905d0, 181.979800d0, 100.466456d0, 355.432999d0, 34.351519d0, 50.077444d0, 314.055005d0, 304.348665d0/)
   real(kind = dp), dimension(8), parameter :: planet1 = (/538101628.68898d0, 210664136.43355d0, 129597742.27580d0, 68905077.59284d0, 10925660.42861d0, 4399609.65932d0, 1542481.19393d0, 786550.32074d0/)

   !
   ! Offest ICRS / equinoxe
   !

   real(kind = dp), parameter :: phid = 0.05542d0/3600d0 !Chapront et al. 2002 gamma to O_icrs in deg ; 0.05542 in arcsec
!   real(kind = 8), parameter :: phid = 0.0d0/3600d0 !Chapront et al. 2002 gamma to O_icrs in deg ; 0.05542 in arcsec



   integer, dimension(4, 218, 3) :: mcent
   real(kind = dp), dimension(218, 3) :: ampli

   integer, dimension(13, 244, 3) :: mplan0
   real(kind = dp), dimension(244, 3) :: phase0, amplip0

   integer, dimension(13, 154, 3) :: mplan1
   real(kind = dp), dimension(154, 3) :: phase1, amplip1

   integer, dimension(13, 25, 3) :: mplan2
   real(kind = dp), dimension(25, 3) :: phase2, amplip2

   real(kind = dp), dimension(13) :: parg
   real(kind = dp), dimension(4) :: darg

   real(kind = dp) :: tt, tt2, tt3, tt4, deltat

   real(kind = dp) :: arg, x(3), y(3)
   real(kind = dp) :: posm(3), dpsi, deps
   ! real(kind=dp)                          :: dist_moon,  epsilon
   real(kind = dp) :: epsilon
   integer :: i, k, l, lmax
   integer :: id, ip

   !
   !  218     Longitude main problem  Order 0 amplitude in arcsec
   !
   !      D   l'  l   F
   DATA (k, (mcent(I, l, 1), I = 1, 4), ampli(l, 1), l = 1, 218) /&
   1, 0, 0, 1, 0, 22639.58578d0, &
   2, 2, 0, -1, 0, 4586.43830d0, &
   3, 2, 0, 0, 0, 2369.91394d0, &
   4, 0, 0, 2, 0, 769.02571d0, &
   5, 0, 1, 0, 0, -666.41710d0, &
   6, 0, 0, 0, 2, -411.59567d0, &
   7, 2, 0, -2, 0, 211.65555d0, &
   8, 2, -1, -1, 0, 205.43582d0, &
   9, 2, 0, 1, 0, 191.95620d0, &
   10, 2, -1, 0, 0, 164.72851d0, &
   11, 0, 1, -1, 0, -147.32129d0, &
   12, 1, 0, 0, 0, -124.98812d0, &
   13, 0, 1, 1, 0, -109.38029d0, &
   14, 2, 0, 0, -2, 55.17705d0, &
   15, 0, 0, 1, 2, -45.09960d0, &
   16, 0, 0, 1, -2, 39.53329d0, &
   17, 4, 0, -1, 0, 38.42983d0, &
   18, 0, 0, 3, 0, 36.12381d0, &
   19, 4, 0, -2, 0, 30.77257d0, &
   20, 2, 1, -1, 0, -28.39708d0, &
   21, 2, 1, 0, 0, -24.35821d0, &
   22, 1, 0, -1, 0, -18.58471d0, &
   23, 1, 1, 0, 0, 17.95446d0, &
   24, 2, -1, 1, 0, 14.53027d0, &
   25, 2, 0, 2, 0, 14.37970d0, &
   26, 4, 0, 0, 0, 13.89906d0, &
   27, 2, 0, -3, 0, 13.19406d0, &
   28, 0, 1, -2, 0, -9.67905d0, &
   29, 2, 0, -1, 2, -9.36586d0, &
   30, 2, -1, -2, 0, 8.60553d0, &
   31, 1, 0, 1, 0, -8.45310d0, &
   32, 2, -2, 0, 0, 8.05016d0, &
   33, 0, 1, 2, 0, -7.63015d0, &
   34, 0, 2, 0, 0, -7.44749d0, &
   35, 2, -2, -1, 0, 7.37119d0, &
   36, 2, 0, 1, -2, -6.38315d0, &
   37, 2, 0, 0, 2, -5.74161d0, &
   38, 4, -1, -1, 0, 4.37401d0, &
   39, 0, 0, 2, 2, -3.99761d0, &
   40, 3, 0, -1, 0, -3.20969d0, &
   41, 2, 1, 1, 0, -2.91454d0, &
   42, 4, -1, -2, 0, 2.73189d0, &
   43, 0, 2, -1, 0, -2.56794d0, &
   44, 2, 2, -1, 0, -2.52120d0, &
   45, 2, 1, -2, 0, 2.48889d0, &
   46, 2, -1, 0, -2, 2.14607d0, &
   47, 4, 0, 1, 0, 1.97773d0, &
   48, 0, 0, 4, 0, 1.93368d0, &
   49, 4, -1, 0, 0, 1.87076d0, &
   50, 1, 0, -2, 0, -1.75297d0, &
   51, 2, 1, 0, -2, -1.43716d0, &
   52, 0, 0, 2, -2, -1.37257d0, &
   53, 1, 1, 1, 0, 1.26182d0, &
   54, 3, 0, -2, 0, -1.22412d0, &
   55, 4, 0, -3, 0, 1.18683d0, &
   56, 2, -1, 2, 0, 1.17700d0, &
   57, 0, 2, 1, 0, -1.16169d0, &
   58, 1, 1, -1, 0, 1.07769d0, &
   59, 2, 0, 3, 0, 1.05950d0, &
   60, 2, 0, 1, 2, -0.99022d0, &
   61, 2, 0, -4, 0, 0.94828d0, &
   62, 2, -2, 1, 0, 0.75168d0, &
   63, 0, 1, -3, 0, -0.66938d0, &
   64, 4, 1, -1, 0, -0.63521d0, &
   65, 1, 0, 2, 0, -0.58399d0, &
   66, 1, 0, 0, -2, -0.58331d0, &
   67, 6, 0, -2, 0, 0.57156d0, &
   68, 2, 0, -2, -2, -0.56064d0, &
   69, 1, -1, 0, 0, -0.55692d0, &
   70, 0, 1, 3, 0, -0.54592d0, &
   71, 2, 0, -2, 2, -0.53571d0, &
   72, 2, -1, -3, 0, 0.47840d0, &
   73, 2, 0, 2, -2, -0.45379d0, &
   74, 2, -1, -1, 2, -0.42622d0, &
   75, 0, 0, 0, 4, 0.42033d0, &
   76, 0, 1, 0, 2, 0.41340d0, &
   77, 3, 0, 0, 0, 0.40423d0, &
   78, 6, 0, -1, 0, 0.39451d0, &
   79, 2, -1, 0, 2, -0.38213d0, &
   80, 2, -1, 1, -2, -0.37451d0, &
   81, 4, 1, -2, 0, -0.35758d0, &
   82, 1, 1, -2, 0, 0.34965d0, &
   83, 2, -3, 0, 0, 0.33979d0, &
   84, 0, 0, 3, 2, -0.32866d0, &
   85, 4, -2, -1, 0, 0.30872d0, &
   86, 0, 1, -1, -2, 0.30155d0, &
   87, 4, 0, -1, -2, 0.30086d0, &
   88, 2, -2, -2, 0, 0.29420d0, &
   89, 6, 0, -3, 0, 0.29255d0, &
   90, 2, 1, 2, 0, -0.29022d0, &
   91, 4, 1, 0, 0, -0.28910d0, &
   92, 4, -1, 1, 0, 0.28250d0, &
   93, 3, 1, -1, 0, 0.27376d0, &
   94, 0, 1, 1, 2, 0.26337d0, &
   95, 1, 0, 0, 2, 0.25429d0, &
   96, 3, 0, 0, -2, -0.25304d0, &
   97, 2, 2, -2, 0, -0.24988d0, &
   98, 2, -3, -1, 0, 0.24694d0, &
   99, 3, -1, -1, 0, -0.23140d0, &
   100, 4, 0, 2, 0, 0.21853d0, &
   101, 4, 0, -1, 2, -0.20134d0, &
   102, 0, 2, -2, 0, -0.19310d0, &
   103, 2, 2, 0, 0, -0.18575d0, &
   104, 2, 0, -1, -2, 0.17903d0, &
   105, 2, 1, -3, 0, 0.17623d0, &
   106, 4, 0, -2, 2, -0.16977d0, &
   107, 4, -2, -2, 0, 0.15780d0, &
   108, 4, -2, 0, 0, 0.15226d0, &
   109, 3, 1, 0, 0, 0.14989d0, &
   110, 1, -1, -1, 0, -0.13636d0, &
   111, 1, 0, -3, 0, -0.12812d0, &
   112, 6, 0, 0, 0, 0.12616d0, &
   113, 2, 0, 2, 2, -0.12386d0, &
   114, 1, -1, 1, 0, -0.12073d0, &
   115, 0, 0, 5, 0, 0.11100d0, &
   116, 0, 3, 0, 0, -0.10135d0, &
   117, 4, -1, -3, 0, 0.09982d0, &
   118, 2, -1, 3, 0, 0.09320d0, &
   119, 1, 1, 2, 0, 0.09205d0, &
   120, 2, 0, -3, -2, -0.09154d0, &
   121, 0, 0, 1, 4, 0.09092d0, &
   122, 6, -1, -2, 0, 0.09033d0, &
   123, 4, 0, 0, 2, -0.08500d0, &
   124, 2, 1, 1, -2, 0.08472d0, &
   125, 3, -1, -2, 0, -0.08311d0, &
   126, 0, 1, 1, -2, -0.08282d0, &
   127, 0, 1, -1, 2, -0.08049d0, &
   128, 0, 0, 1, -4, -0.08019d0, &
   129, 2, 0, 4, 0, 0.07765d0, &
   130, 2, 0, 0, -4, -0.07518d0, &
   131, 0, 1, 0, -2, 0.07501d0, &
   132, 2, -1, 1, 2, -0.07373d0, &
   133, 6, -1, -1, 0, 0.07142d0, &
   134, 2, 0, -5, 0, 0.06850d0, &
   135, 2, 1, -1, 2, 0.06742d0, &
   136, 4, 0, 1, -2, -0.06601d0, &
   137, 2, 1, 0, 2, 0.06541d0, &
   138, 0, 2, 2, 0, -0.06513d0, &
   139, 3, -1, 0, 0, 0.06507d0, &
   140, 2, -2, 2, 0, 0.06439d0, &
   141, 2, -2, 0, -2, 0.06313d0, &
   142, 2, -1, -1, -2, -0.06103d0, &
   143, 5, 0, -2, 0, -0.05725d0, &
   144, 0, 0, 3, -2, -0.05684d0, &
   145, 0, 1, -2, -2, 0.05165d0, &
   146, 0, 3, -1, 0, -0.05141d0, &
   147, 4, 1, 1, 0, -0.05070d0, &
   148, 0, 1, -4, 0, -0.04702d0, &
   149, 1, 0, 1, 2, 0.04450d0, &
   150, 3, 0, -3, 0, -0.04442d0, &
   151, 0, 1, 2, 2, 0.04338d0, &
   152, 1, -2, 0, 0, 0.04304d0, &
   153, 3, 1, -2, 0, -0.04189d0, &
   154, 1, 0, 3, 0, -0.04074d0, &
   155, 1, 0, 1, -2, -0.04012d0, &
   156, 1, 2, 0, 0, -0.03968d0, &
   157, 0, 1, 4, 0, -0.03947d0, &
   158, 6, -1, -3, 0, 0.03900d0, &
   159, 1, 1, 0, 2, -0.03587d0, &
   160, 4, 2, -2, 0, -0.03514d0, &
   161, 2, 0, 3, -2, -0.03336d0, &
   162, 2, -3, 1, 0, 0.03300d0, &
   163, 4, -1, 2, 0, 0.03274d0, &
   164, 3, 0, -1, -2, -0.02979d0, &
   165, 2, -1, -4, 0, 0.02949d0, &
   166, 2, -1, 2, -2, -0.02887d0, &
   167, 2, -1, -2, -2, -0.02804d0, &
   168, 4, 1, -3, 0, 0.02682d0, &
   169, 0, 1, 2, -2, 0.02677d0, &
   170, 2, 1, 3, 0, -0.02676d0, &
   171, 0, 0, 4, 2, -0.02602d0, &
   172, 6, -1, 0, 0, 0.02510d0, &
   173, 0, 1, -2, 2, 0.02429d0, &
   174, 4, -2, 1, 0, 0.02411d0, &
   175, 4, 0, 0, -2, -0.02391d0, &
   176, 1, 0, -1, -2, -0.02379d0, &
   177, 2, 2, 0, -2, -0.02349d0, &
   178, 1, 1, -3, 0, 0.02296d0, &
   179, 4, -1, -1, -2, 0.02289d0, &
   180, 6, 0, 1, 0, 0.02285d0, &
   181, 4, -1, -1, 2, -0.02273d0, &
   182, 3, 1, 1, 0, 0.02244d0, &
   183, 4, 2, -1, 0, -0.02171d0, &
   184, 2, -1, -2, 2, -0.02157d0, &
   185, 4, 0, 3, 0, 0.02149d0, &
   186, 2, 0, -1, 4, 0.01993d0, &
   187, 3, -1, 0, -2, -0.01948d0, &
   188, 4, 1, -1, -2, -0.01875d0, &
   189, 2, 1, -4, 0, 0.01819d0, &
   190, 2, -2, 0, 2, -0.01816d0, &
   191, 0, 3, 1, 0, -0.01796d0, &
   192, 4, 0, 1, 2, -0.01781d0, &
   193, 4, -3, -1, 0, 0.01741d0, &
   194, 5, 0, -3, 0, -0.01686d0, &
   195, 2, -2, 1, -2, -0.01644d0, &
   196, 2, 1, 1, 2, 0.01605d0, &
   197, 1, 0, -1, 2, 0.01598d0, &
   198, 2, -2, -3, 0, 0.01544d0, &
   199, 2, -2, -1, 2, -0.01541d0, &
   200, 4, -1, -2, 2, -0.01533d0, &
   201, 0, 2, -3, 0, -0.01514d0, &
   202, 1, -1, 2, 0, -0.01483d0, &
   203, 6, 0, -4, 0, 0.01376d0, &
   204, 2, 0, 0, 4, 0.01372d0, &
   205, 5, 0, -1, 0, -0.01350d0, &
   206, 2, 2, 1, 0, -0.01343d0, &
   207, 2, 0, 3, 2, -0.01332d0, &
   208, 2, -4, 0, 0, 0.01331d0, &
   209, 0, 0, 2, 4, 0.01297d0, &
   210, 6, 1, -2, 0, -0.01282d0, &
   211, 1, -1, 0, -2, -0.01281d0, &
   212, 3, 0, -1, 2, 0.01215d0, &
   213, 3, -2, -1, 0, -0.01182d0, &
   214, 4, -1, 0, 2, -0.01114d0, &
   215, 2, 0, -4, -2, -0.01077d0, &
   216, 6, 1, -1, 0, -0.01064d0, &
   217, 3, 0, 1, -2, -0.01062d0, &
   218, 2, -1, 2, 2, -0.01007d0 /
   !
   ! 188     Latitude  main problem  Order 0
   !
   !      D   l'  l   F
   DATA (k, (mcent(I, l, 2), I = 1, 4), ampli(l, 2), l = 1, 188) /&
   1, 0, 0, 0, 1, 18461.23868d0, &
   2, 0, 0, 1, 1, 1010.16707d0, &
   3, 0, 0, 1, -1, 999.69358d0, &
   4, 2, 0, 0, -1, 623.65243d0, &
   5, 2, 0, -1, 1, 199.48374d0, &
   6, 2, 0, -1, -1, 166.57410d0, &
   7, 2, 0, 0, 1, 117.26069d0, &
   8, 0, 0, 2, 1, 61.91195d0, &
   9, 2, 0, 1, -1, 33.35720d0, &
   10, 0, 0, 2, -1, 31.75967d0, &
   11, 2, -1, 0, -1, 29.57658d0, &
   12, 2, 0, -2, -1, 15.56626d0, &
   13, 2, 0, 1, 1, 15.12155d0, &
   14, 2, 1, 0, -1, -12.09414d0, &
   15, 2, -1, -1, 1, 8.86814d0, &
   16, 2, -1, 0, 1, 7.95855d0, &
   17, 2, -1, -1, -1, 7.43455d0, &
   18, 0, 1, -1, -1, -6.73143d0, &
   19, 4, 0, -1, -1, 6.57957d0, &
   20, 0, 1, 0, 1, -6.46007d0, &
   21, 0, 0, 0, 3, -6.29648d0, &
   22, 0, 1, -1, 1, -5.63235d0, &
   23, 1, 0, 0, 1, -5.36840d0, &
   24, 0, 1, 1, 1, -5.31127d0, &
   25, 0, 1, 1, -1, -5.07591d0, &
   26, 0, 1, 0, -1, -4.83961d0, &
   27, 1, 0, 0, -1, -4.80574d0, &
   28, 0, 0, 3, 1, 3.98405d0, &
   29, 4, 0, 0, -1, 3.67446d0, &
   30, 4, 0, -1, 1, 2.99848d0, &
   31, 0, 0, 1, -3, 2.79864d0, &
   32, 4, 0, -2, 1, 2.41388d0, &
   33, 2, 0, 0, -3, 2.18631d0, &
   34, 2, 0, 2, -1, 2.14617d0, &
   35, 2, -1, 1, -1, 1.76598d0, &
   36, 2, 0, -2, 1, -1.62442d0, &
   37, 0, 0, 3, -1, 1.58130d0, &
   38, 2, 0, 2, 1, 1.51975d0, &
   39, 2, 0, -3, -1, 1.51563d0, &
   40, 2, 1, -1, 1, -1.31782d0, &
   41, 2, 1, 0, 1, -1.26427d0, &
   42, 4, 0, 0, 1, 1.19187d0, &
   43, 2, -1, 1, 1, 1.13461d0, &
   44, 2, -2, 0, -1, 1.08578d0, &
   45, 0, 0, 1, 3, -1.01938d0, &
   46, 2, 1, 1, -1, -0.82271d0, &
   47, 1, 1, 0, -1, 0.80422d0, &
   48, 1, 1, 0, 1, 0.80259d0, &
   49, 0, 1, -2, -1, -0.79319d0, &
   50, 2, 1, -1, -1, -0.79101d0, &
   51, 1, 0, 1, 1, -0.66741d0, &
   52, 2, -1, -2, -1, 0.65022d0, &
   53, 0, 1, 2, 1, -0.63881d0, &
   54, 4, 0, -2, -1, 0.63371d0, &
   55, 4, -1, -1, -1, 0.59577d0, &
   56, 1, 0, 1, -1, -0.58893d0, &
   57, 4, 0, 1, -1, 0.47338d0, &
   58, 1, 0, -1, -1, -0.42989d0, &
   59, 4, -1, 0, -1, 0.41494d0, &
   60, 2, -2, 0, 1, 0.38350d0, &
   61, 3, 0, 0, -1, -0.35183d0, &
   62, 4, -1, -1, 1, 0.33881d0, &
   63, 2, 0, -1, -3, 0.32906d0, &
   64, 2, -2, -1, 1, 0.31471d0, &
   65, 0, 1, 2, -1, -0.31291d0, &
   66, 3, 0, -1, -1, -0.30517d0, &
   67, 0, 1, -2, 1, -0.30128d0, &
   68, 2, 0, 1, -3, -0.29115d0, &
   69, 2, -2, -1, -1, 0.26863d0, &
   70, 0, 0, 4, 1, 0.26325d0, &
   71, 2, 0, -3, 1, 0.25408d0, &
   72, 2, 0, -1, 3, -0.24483d0, &
   73, 2, 1, 1, 1, -0.23701d0, &
   74, 4, -1, -2, 1, 0.21375d0, &
   75, 4, 0, 1, 1, 0.21259d0, &
   76, 3, 0, -1, 1, -0.20593d0, &
   77, 4, 1, -1, -1, -0.17190d0, &
   78, 4, -1, 0, 1, 0.15790d0, &
   79, 2, 0, 3, -1, 0.14642d0, &
   80, 2, 0, 0, 3, -0.14453d0, &
   81, 1, 0, -1, 1, 0.13928d0, &
   82, 2, 0, 3, 1, 0.13795d0, &
   83, 2, 2, 0, -1, -0.13414d0, &
   84, 2, 0, -4, -1, 0.13381d0, &
   85, 0, 0, 2, -3, -0.13035d0, &
   86, 2, -1, 2, -1, 0.12896d0, &
   87, 2, -1, 2, 1, 0.12386d0, &
   88, 0, 0, 2, 3, -0.11787d0, &
   89, 0, 2, -1, -1, -0.11334d0, &
   90, 2, 2, -1, 1, -0.11329d0, &
   91, 4, 1, 0, -1, -0.11307d0, &
   92, 1, 0, -2, -1, -0.10964d0, &
   93, 2, 2, -1, -1, -0.10534d0, &
   94, 1, 1, 1, 1, 0.10176d0, &
   95, 0, 2, -1, 1, -0.09510d0, &
   96, 6, 0, -1, -1, 0.09403d0, &
   97, 0, 0, 4, -1, 0.09157d0, &
   98, 2, -1, 0, -3, 0.08814d0, &
   99, 6, 0, -2, -1, 0.08096d0, &
   100, 2, 1, -2, -1, 0.07913d0, &
   101, 1, 0, -2, 1, -0.07846d0, &
   102, 0, 1, -3, -1, -0.07479d0, &
   103, 2, -2, 1, -1, 0.06914d0, &
   104, 2, 0, -2, 3, -0.06561d0, &
   105, 1, 0, 2, 1, -0.06383d0, &
   106, 2, 1, 2, -1, -0.06283d0, &
   107, 4, 0, 0, -3, 0.06257d0, &
   108, 2, -1, -2, 1, -0.06208d0, &
   109, 0, 2, 1, -1, -0.06186d0, &
   110, 0, 1, 3, 1, -0.06176d0, &
   111, 6, 0, -2, 1, 0.05963d0, &
   112, 2, -2, 1, 1, 0.05848d0, &
   113, 0, 2, 0, 1, -0.05729d0, &
   114, 4, -1, 1, -1, 0.05686d0, &
   115, 1, 1, -1, 1, -0.05590d0, &
   116, 0, 2, 1, 1, -0.05504d0, &
   117, 2, -1, -3, -1, 0.05502d0, &
   118, 2, 1, 0, -3, -0.05457d0, &
   119, 2, 1, -2, 1, 0.05429d0, &
   120, 4, -1, -2, -1, 0.05251d0, &
   121, 4, 1, -1, 1, -0.05097d0, &
   122, 3, 0, -2, 1, -0.04852d0, &
   123, 4, 0, 2, -1, 0.04834d0, &
   124, 6, 0, -1, 1, 0.04217d0, &
   125, 3, 0, -2, -1, -0.03941d0, &
   126, 6, 0, 0, -1, 0.03674d0, &
   127, 2, -3, 0, -1, 0.03647d0, &
   128, 1, 0, 2, -1, -0.03636d0, &
   129, 3, 0, 1, -1, -0.03611d0, &
   130, 1, 1, 1, -1, 0.03465d0, &
   131, 4, -2, -1, -1, 0.03462d0, &
   132, 3, 1, 0, -1, 0.03436d0, &
   133, 1, 0, 0, -3, -0.03226d0, &
   134, 2, 1, 2, 1, -0.03142d0, &
   135, 6, 0, -3, 1, 0.03118d0, &
   136, 2, 0, 1, 3, -0.03038d0, &
   137, 4, -1, 1, 1, 0.03009d0, &
   138, 4, 1, -2, 1, -0.02957d0, &
   139, 4, -2, 0, -1, 0.02899d0, &
   140, 3, 0, 0, 1, -0.02840d0, &
   141, 4, 0, 2, 1, 0.02828d0, &
   142, 3, -1, 0, -1, -0.02572d0, &
   143, 4, 1, 0, 1, -0.02549d0, &
   144, 2, 0, -4, 1, 0.02496d0, &
   145, 0, 1, 3, -1, -0.02419d0, &
   146, 4, -2, -1, 1, 0.02380d0, &
   147, 0, 1, -3, 1, -0.02365d0, &
   148, 2, -2, -2, -1, 0.02285d0, &
   149, 4, 0, -3, 1, 0.02174d0, &
   150, 3, -1, -1, -1, -0.02104d0, &
   151, 3, 1, -1, 1, 0.02083d0, &
   152, 2, 0, -2, -3, 0.02045d0, &
   153, 1, -1, 1, -1, -0.02012d0, &
   154, 1, -1, 0, 1, -0.01829d0, &
   155, 2, 1, -3, -1, 0.01818d0, &
   156, 0, 2, 0, -1, -0.01801d0, &
   157, 0, 0, 5, 1, 0.01768d0, &
   158, 4, 1, 1, -1, -0.01692d0, &
   159, 1, 1, -2, 1, 0.01680d0, &
   160, 2, -1, 1, -3, -0.01669d0, &
   161, 2, -3, 0, 1, 0.01603d0, &
   162, 1, 1, -2, -1, 0.01597d0, &
   163, 0, 2, -2, -1, -0.01571d0, &
   164, 6, -1, -1, -1, 0.01486d0, &
   165, 2, 2, 0, 1, -0.01482d0, &
   166, 6, 0, 0, 1, 0.01465d0, &
   167, 3, -1, -1, 1, -0.01356d0, &
   168, 3, 1, 0, 1, 0.01351d0, &
   169, 1, -1, 0, -1, -0.01346d0, &
   170, 3, 1, -1, -1, 0.01321d0, &
   171, 4, -2, 0, 1, 0.01270d0, &
   172, 2, 2, -2, -1, -0.01262d0, &
   173, 1, 0, -3, -1, -0.01255d0, &
   174, 4, -2, -2, 1, 0.01230d0, &
   175, 2, -1, 3, 1, 0.01211d0, &
   176, 2, 0, 4, 1, 0.01186d0, &
   177, 0, 0, 3, 3, -0.01181d0, &
   178, 2, -1, -1, 3, -0.01177d0, &
   179, 0, 1, 0, 3, 0.01157d0, &
   180, 2, 0, -5, -1, 0.01127d0, &
   181, 6, -1, -2, -1, 0.01091d0, &
   182, 5, 0, -1, -1, -0.01049d0, &
   183, 2, -3, -1, 1, 0.01042d0, &
   184, 2, -1, -1, -3, 0.01034d0, &
   185, 1, -1, -1, -1, 0.01031d0, &
   186, 2, 0, 4, -1, 0.01027d0, &
   187, 1, 1, 2, 1, 0.01016d0, &
   188, 3, 0, 0, -3, -0.01009d0 /
   !
   ! 155     Distance  main problem  order 0
   !
   !      D   l'  l   F
   DATA (k, (mcent(I, l, 3), I = 1, 4), ampli(l, 3), l = 1, 155) /&
   1, 0, 0, 0, 0, 385000.52899d0, &
   2, 0, 0, 1, 0, -20905.35504d0, &
   3, 2, 0, -1, 0, -3699.11092d0, &
   4, 2, 0, 0, 0, -2955.96756d0, &
   5, 0, 0, 2, 0, -569.92512d0, &
   6, 2, 0, -2, 0, 246.15848d0, &
   7, 2, -1, 0, 0, -204.58598d0, &
   8, 2, 0, 1, 0, -170.73308d0, &
   9, 2, -1, -1, 0, -152.13771d0, &
   10, 0, 1, -1, 0, -129.62014d0, &
   11, 1, 0, 0, 0, 108.74270d0, &
   12, 0, 1, 1, 0, 104.75523d0, &
   13, 0, 0, 1, -2, 79.66056d0, &
   14, 0, 1, 0, 0, 48.88830d0, &
   15, 4, 0, -1, 0, -34.78252d0, &
   16, 2, 1, 0, 0, 30.82384d0, &
   17, 2, 1, -1, 0, 24.20848d0, &
   18, 0, 0, 3, 0, -23.21043d0, &
   19, 4, 0, -2, 0, -21.63634d0, &
   20, 1, 1, 0, 0, -16.67471d0, &
   21, 2, 0, -3, 0, 14.40269d0, &
   22, 2, -1, 1, 0, -12.83140d0, &
   23, 4, 0, 0, 0, -11.64995d0, &
   24, 2, 0, 2, 0, -10.44476d0, &
   25, 2, 0, 0, -2, 10.32111d0, &
   26, 2, -1, -2, 0, 10.05620d0, &
   27, 2, -2, 0, 0, -9.88445d0, &
   28, 2, 0, -1, -2, 8.75156d0, &
   29, 1, 0, -1, 0, -8.37911d0, &
   30, 0, 1, -2, 0, -7.00269d0, &
   31, 1, 0, 1, 0, 6.32200d0, &
   32, 0, 1, 2, 0, 5.75085d0, &
   33, 2, -2, -1, 0, -4.95013d0, &
   34, 0, 0, 2, -2, -4.42118d0, &
   35, 2, 0, 1, -2, 4.13111d0, &
   36, 4, -1, -1, 0, -3.95798d0, &
   37, 3, 0, -1, 0, 3.25824d0, &
   38, 0, 0, 0, 2, -3.14830d0, &
   39, 2, 1, 1, 0, 2.61641d0, &
   40, 2, 2, -1, 0, 2.35363d0, &
   41, 0, 2, -1, 0, -2.11713d0, &
   42, 4, -1, -2, 0, -1.89704d0, &
   43, 1, 0, -2, 0, -1.73853d0, &
   44, 4, -1, 0, 0, -1.57139d0, &
   45, 4, 0, 1, 0, -1.42255d0, &
   46, 3, 0, 0, 0, -1.41893d0, &
   47, 0, 2, 1, 0, 1.16553d0, &
   48, 0, 0, 4, 0, -1.11694d0, &
   49, 0, 2, 0, 0, 1.06567d0, &
   50, 1, 1, 1, 0, -0.93332d0, &
   51, 3, 0, -2, 0, 0.86243d0, &
   52, 1, 1, -1, 0, 0.85124d0, &
   53, 2, -1, 2, 0, -0.84880d0, &
   54, 1, 0, 0, -2, -0.79563d0, &
   55, 2, 0, -4, 0, 0.77854d0, &
   56, 2, 0, -2, 2, 0.77404d0, &
   57, 2, 0, 3, 0, -0.66968d0, &
   58, 2, -2, 1, 0, -0.65753d0, &
   59, 2, -1, 0, -2, 0.65706d0, &
   60, 2, 0, -1, 2, 0.59632d0, &
   61, 4, 1, -1, 0, 0.57879d0, &
   62, 4, 0, -3, 0, -0.51423d0, &
   63, 4, 0, 0, -2, -0.50792d0, &
   64, 1, -1, 0, 0, 0.49755d0, &
   65, 2, -1, -3, 0, 0.49504d0, &
   66, 2, 0, -2, -2, 0.47262d0, &
   67, 6, 0, -2, 0, -0.42250d0, &
   68, 0, 1, -3, 0, -0.42241d0, &
   69, 2, -3, 0, 0, -0.41071d0, &
   70, 1, 0, 2, 0, 0.37852d0, &
   71, 0, 1, 3, 0, 0.35508d0, &
   72, 2, -2, -2, 0, 0.34302d0, &
   73, 0, 1, -1, 2, 0.33463d0, &
   74, 1, 1, -2, 0, 0.33225d0, &
   75, 2, -1, -1, -2, 0.32334d0, &
   76, 4, 0, -1, -2, -0.32176d0, &
   77, 6, 0, -1, 0, -0.28663d0, &
   78, 2, 0, 2, -2, 0.28399d0, &
   79, 4, -2, -1, 0, -0.27904d0, &
   80, 3, -1, -1, 0, 0.25560d0, &
   81, 0, 1, 1, -2, -0.24810d0, &
   82, 4, 1, 0, 0, 0.24452d0, &
   83, 4, 1, -2, 0, 0.23695d0, &
   84, 3, 1, -1, 0, -0.21258d0, &
   85, 2, 1, 2, 0, 0.21251d0, &
   86, 2, -1, 1, -2, 0.20941d0, &
   87, 4, -1, 1, 0, -0.20285d0, &
   88, 3, 0, 0, -2, 0.20099d0, &
   89, 0, 1, 0, -2, -0.18567d0, &
   90, 6, 0, -3, 0, -0.18316d0, &
   91, 2, 1, -3, 0, 0.16857d0, &
   92, 0, 1, 0, 2, -0.15802d0, &
   93, 3, -1, 0, 0, -0.15707d0, &
   94, 2, -3, -1, 0, -0.14806d0, &
   95, 2, 2, 0, 0, 0.14763d0, &
   96, 2, 1, -2, 0, 0.14368d0, &
   97, 4, 0, 2, 0, -0.13922d0, &
   98, 0, 2, -2, 0, -0.13617d0, &
   99, 2, 1, 0, -2, -0.13571d0, &
   100, 4, -2, 0, 0, -0.12805d0, &
   101, 1, -1, -1, 0, 0.11411d0, &
   102, 1, -1, 1, 0, 0.10998d0, &
   103, 2, 2, -2, 0, -0.10887d0, &
   104, 4, -2, -2, 0, -0.10833d0, &
   105, 3, 1, 0, 0, -0.10766d0, &
   106, 0, 0, 1, 2, -0.10326d0, &
   107, 1, 0, -3, 0, -0.09938d0, &
   108, 6, 0, 0, 0, -0.08587d0, &
   109, 4, 0, -2, -2, -0.07982d0, &
   110, 6, -1, -2, 0, -0.06678d0, &
   111, 3, 0, 1, 0, -0.06545d0, &
   112, 1, 0, 1, -2, 0.06055d0, &
   113, 1, 1, 2, 0, -0.05904d0, &
   114, 0, 0, 5, 0, -0.05888d0, &
   115, 2, -1, 3, 0, -0.05850d0, &
   116, 4, -1, 0, -2, -0.05789d0, &
   117, 2, 1, 1, -2, -0.05527d0, &
   118, 3, -1, -2, 0, 0.05293d0, &
   119, 6, -1, -1, 0, -0.05191d0, &
   120, 0, 2, 2, 0, 0.05072d0, &
   121, 0, 1, -2, 2, -0.05020d0, &
   122, 3, 0, -3, 0, -0.04843d0, &
   123, 2, 0, -5, 0, 0.04740d0, &
   124, 2, 1, -1, -2, -0.04736d0, &
   125, 2, -2, 2, 0, -0.04608d0, &
   126, 5, 0, -2, 0, 0.04591d0, &
   127, 2, 0, 4, 0, -0.04422d0, &
   128, 4, -1, -3, 0, -0.04316d0, &
   129, 1, 0, -1, -2, -0.04232d0, &
   130, 0, 3, -1, 0, -0.03894d0, &
   131, 3, 1, -2, 0, 0.03810d0, &
   132, 2, -1, -1, 2, 0.03734d0, &
   133, 1, 2, 0, 0, 0.03729d0, &
   134, 4, 1, 1, 0, 0.03682d0, &
   135, 1, 1, 0, -2, 0.03379d0, &
   136, 0, 1, 2, -2, 0.03265d0, &
   137, 2, 0, 0, 2, 0.03143d0, &
   138, 2, -1, -2, 2, 0.03024d0, &
   139, 1, -2, 0, 0, -0.02948d0, &
   140, 4, 0, -4, 0, -0.02939d0, &
   141, 2, 0, -3, -2, 0.02910d0, &
   142, 2, -3, 1, 0, -0.02855d0, &
   143, 2, -2, 0, -2, 0.02839d0, &
   144, 4, -1, -1, -2, -0.02698d0, &
   145, 0, 1, -4, 0, -0.02674d0, &
   146, 4, 2, -2, 0, 0.02658d0, &
   147, 1, 0, -1, 2, -0.02471d0, &
   148, 6, -1, -3, 0, -0.02436d0, &
   149, 4, 1, -3, 0, -0.02399d0, &
   150, 1, 0, 3, 0, 0.02368d0, &
   151, 2, -1, -4, 0, 0.02334d0, &
   152, 0, 1, 4, 0, 0.02304d0, &
   153, 0, 3, 0, 0, 0.02127d0, &
   154, 4, -1, 2, 0, -0.02079d0, &
   155, 2, 0, -3, 2, -0.02008d0 /
   !
   ! 244     Longitude Perturbations     order 0
   !
   DATA (k, (mplan0(I, l, 1), I = 1, 13), phase0(l, 1), amplip0(l, 1), l = 1, 244) /&
   1, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 26.54261d0, 14.24883d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0.00094d0, 7.06304d0, &
   3, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -1, 0, 180.11977d0, 1.14307d0, &
   4, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 285.98707d0, 0.90114d0, &
   5, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 180.00988d0, 0.82155d0, &
   6, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 26.54324d0, 0.78811d0, &
   7, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26.54560d0, 0.73930d0, &
   8, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 179.98144d0, 0.64371d0, &
   9, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1.22890d0, 0.63880d0, &
   10, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 333.30551d0, 0.56341d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -1, 0.00127d0, 0.49331d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0.00127d0, 0.49141d0, &
   13, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 10.07001d0, 0.44532d0, &
   14, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0.00071d0, 0.36061d0, &
   15, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 269.95393d0, 0.34355d0, &
   16, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 318.13776d0, 0.32455d0, &
   17, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.20448d0, 0.30155d0, &
   18, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 95.13523d0, 0.28938d0, &
   19, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -2, 0, 10.03835d0, 0.28281d0, &
   20, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -2, 0, 0.08642d0, 0.24515d0, &
   21, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 235.75013d0, 0.22653d0, &
   22, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, 1, 0, 1.74333d0, 0.21118d0, &
   23, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 359.76487d0, 0.19443d0, &
   24, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 180.43608d0, 0.18457d0, &
   25, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, -1, 0, 359.98450d0, 0.18256d0, &
   26, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 270.99316d0, 0.17511d0, &
   27, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, -1, 0, 284.98776d0, 0.17083d0, &
   28, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 1, 0, 284.98767d0, 0.17082d0, &
   29, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 26.53985d0, 0.16710d0, &
   30, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, 0, 0, 1.80280d0, 0.16497d0, &
   31, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 26.54223d0, 0.16440d0, &
   32, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 301.00599d0, 0.16425d0, &
   33, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 0, 0, 347.76653d0, 0.16376d0, &
   34, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 26.54301d0, 0.16302d0, &
   35, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 0.31177d0, 0.16090d0, &
   36, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 26.54227d0, 0.15711d0, &
   37, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 180.00888d0, 0.15350d0, &
   38, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 1.23582d0, 0.14347d0, &
   39, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, -1, 0, 287.11388d0, 0.13971d0, &
   40, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 1, 0, 287.11216d0, 0.13960d0, &
   41, 0, 2, -2, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 359.99762d0, 0.13593d0, &
   42, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0.36659d0, 0.13432d0, &
   43, 0, 1, -1, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 179.98436d0, 0.13122d0, &
   44, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 180.01511d0, 0.12722d0, &
   45, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0.12046d0, 0.12537d0, &
   46, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 1.11168d0, 0.10993d0, &
   47, 0, 20, -21, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 271.80227d0, 0.10651d0, &
   48, 0, 26, -29, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 75.00074d0, 0.10489d0, &
   49, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 270.34292d0, 0.10386d0, &
   50, 0, 1, -1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 179.98152d0, 0.09922d0, &
   51, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, -1, 0.00090d0, 0.09642d0, &
   52, 0, 0, 2, 0, -3, 0, 0, 0, 0, 0, 0, -1, 0, 190.17970d0, 0.09518d0, &
   53, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 300.13913d0, 0.09178d0, &
   54, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, 0, 0, 286.69269d0, 0.09009d0, &
   55, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 0, 0, 286.68957d0, 0.09002d0, &
   56, 0, 0, 2, 0, -2, 0, 0, 0, 0, -2, 0, 1, 0, 180.86110d0, 0.08911d0, &
   57, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 286.69628d0, 0.08634d0, &
   58, 0, 0, 2, 0, -2, 0, 0, 0, 0, -2, 0, 0, 0, 180.95918d0, 0.08408d0, &
   59, 0, 3, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 344.52995d0, 0.08277d0, &
   60, 0, 4, -4, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 180.00979d0, 0.08235d0, &
   61, 0, 6, -8, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 167.22262d0, 0.07667d0, &
   62, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0.04167d0, 0.07199d0, &
   63, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, 0, 0, 180.13179d0, 0.06863d0, &
   64, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 1, 0, 346.48972d0, 0.06826d0, &
   65, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -1, 0, 346.48704d0, 0.06826d0, &
   66, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 180.70913d0, 0.06711d0, &
   67, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, -1, 0.00100d0, 0.06569d0, &
   68, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 270.28604d0, 0.06514d0, &
   69, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, -1, -1, 0.00042d0, 0.06456d0, &
   70, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 313.47920d0, 0.06406d0, &
   71, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, -1, 163.72199d0, 0.06369d0, &
   72, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 313.30347d0, 0.06323d0, &
   73, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1, 95.13453d0, 0.06245d0, &
   74, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 271.56625d0, 0.06188d0, &
   75, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 95.13459d0, 0.06174d0, &
   76, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 95.13225d0, 0.06154d0, &
   77, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 268.69648d0, 0.06030d0, &
   78, 0, 6, -8, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 167.12686d0, 0.05946d0, &
   79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0.07666d0, 0.05932d0, &
   80, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0, 1, 0, 358.86873d0, 0.05904d0, &
   81, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 0, 1, 0, 359.80946d0, 0.05756d0, &
   82, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0.39852d0, 0.05635d0, &
   83, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 268.56875d0, 0.05540d0, &
   84, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 26.54324d0, 0.05354d0, &
   85, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, -1, 0, 356.13718d0, 0.05274d0, &
   86, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 272.29577d0, 0.05272d0, &
   87, 0, 0, 2, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 190.43600d0, 0.05167d0, &
   88, 0, 0, 1, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 0.57920d0, 0.05108d0, &
   89, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0.15681d0, 0.05071d0, &
   90, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1, 0.00051d0, 0.05036d0, &
   91, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 26.54560d0, 0.05022d0, &
   92, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 1, -1, 0.00029d0, 0.04962d0, &
   93, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 1, 0.00076d0, 0.04746d0, &
   94, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 330.08714d0, 0.04667d0, &
   95, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 270.60642d0, 0.04662d0, &
   96, 0, 4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 180.07081d0, 0.04626d0, &
   97, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 358.57975d0, 0.04612d0, &
   98, 0, 0, 1, 0, -2, 0, 0, 0, 0, -2, 0, 1, 0, 305.99416d0, 0.04570d0, &
   99, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 103.71567d0, 0.04478d0, &
   100, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 179.09975d0, 0.04340d0, &
   101, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2, 0, -1, 0, 151.73240d0, 0.04282d0, &
   102, 0, 20, -20, 0, 0, 0, 0, 0, 0, -1, 0, 1, -1, 194.02829d0, 0.04280d0, &
   103, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 271.38094d0, 0.04214d0, &
   104, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 359.51637d0, 0.04213d0, &
   105, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 294.98414d0, 0.04046d0, &
   106, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 293.67110d0, 0.04040d0, &
   107, 0, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 330.28695d0, 0.04014d0, &
   108, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 271.17391d0, 0.03945d0, &
   109, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0.43296d0, 0.03897d0, &
   110, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 233.44625d0, 0.03878d0, &
   111, 0, 8, -13, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 238.00817d0, 0.03878d0, &
   112, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 233.44527d0, 0.03874d0, &
   113, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 237.96935d0, 0.03850d0, &
   114, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0.00070d0, 0.03838d0, &
   115, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, -1, 0, 301.00167d0, 0.03780d0, &
   116, 0, 5, -5, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 180.07299d0, 0.03759d0, &
   117, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 270.55741d0, 0.03668d0, &
   118, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 180.00000d0, 0.03638d0, &
   119, 0, 2, -2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 359.98667d0, 0.03553d0, &
   120, 0, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, -1, 0, 321.96125d0, 0.03549d0, &
   121, 0, 12, -8, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 228.50839d0, 0.03489d0, &
   122, 0, 0, 1, 0, -2, 0, 0, 0, 0, -2, 0, 0, 0, 306.02804d0, 0.03474d0, &
   123, 0, 2, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 271.01820d0, 0.03451d0, &
   124, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, -1, 0.00126d0, 0.03402d0, &
   125, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 1, 0, 180.40927d0, 0.03337d0, &
   126, 0, 0, 6, -8, 0, 0, 0, 0, 0, 2, 0, -1, 0, 302.34346d0, 0.03317d0, &
   127, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 0, 0, 0, 355.13454d0, 0.03309d0, &
   128, 0, 5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 180.02851d0, 0.03305d0, &
   129, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, -1, 0.00128d0, 0.03279d0, &
   130, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0.07117d0, 0.03259d0, &
   131, 0, 3, -7, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 148.23103d0, 0.03245d0, &
   132, 0, 21, -21, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 179.94393d0, 0.03238d0, &
   133, 0, 0, 1, -2, 0, 0, 0, 0, 0, -2, 0, 1, 0, 320.91988d0, 0.03194d0, &
   134, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, 0, 0, 10.08580d0, 0.03174d0, &
   135, 0, 0, 1, 0, -1, 0, 0, 0, 0, 2, 0, 0, 0, 0.43396d0, 0.03111d0, &
   136, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 333.03408d0, 0.03111d0, &
   137, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 333.10628d0, 0.03091d0, &
   138, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 179.94277d0, 0.03030d0, &
   139, 0, 4, -5, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 91.39917d0, 0.02952d0, &
   140, 3, 0, -1, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 96.56536d0, 0.02923d0, &
   141, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 271.08666d0, 0.02911d0, &
   142, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 1, 0, 301.49768d0, 0.02889d0, &
   143, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 346.26872d0, 0.02848d0, &
   144, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, 206.54280d0, 0.02810d0, &
   145, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, 206.54280d0, 0.02810d0, &
   146, 0, 0, 8, -15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 155.30914d0, 0.02611d0, &
   147, 0, 4, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 343.21324d0, 0.02587d0, &
   148, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, -1, 75.47182d0, 0.02540d0, &
   149, 0, 15, -13, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 26.32315d0, 0.02514d0, &
   150, 0, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 320.19758d0, 0.02513d0, &
   151, 0, 8, -13, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 237.41676d0, 0.02418d0, &
   152, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 237.22280d0, 0.02393d0, &
   153, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 0, 0.69736d0, 0.02338d0, &
   154, 0, 0, 1, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 304.60699d0, 0.02314d0, &
   155, 0, 0, 1, -2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 318.22943d0, 0.02301d0, &
   156, 0, 3, -5, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 346.95020d0, 0.02299d0, &
   157, 0, 5, -8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 70.23554d0, 0.02296d0, &
   158, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 271.88133d0, 0.02267d0, &
   159, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 180.02165d0, 0.02263d0, &
   160, 0, 1, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 180.01464d0, 0.02201d0, &
   161, 0, 2, -3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 269.98647d0, 0.02172d0, &
   162, 0, 0, 3, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 173.71802d0, 0.02166d0, &
   163, 0, 6, -6, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 180.04574d0, 0.02116d0, &
   164, 0, 0, 6, -8, 0, 0, 0, 0, 0, 2, 0, -2, 0, 302.35646d0, 0.02023d0, &
   165, 0, 0, 3, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 300.80387d0, 0.02009d0, &
   166, 0, 0, 4, -4, 0, 0, 0, 0, 0, 2, 0, -1, 0, 181.91276d0, 0.02008d0, &
   167, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, 26.54274d0, 0.01999d0, &
   168, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 26.54486d0, 0.01946d0, &
   169, 0, 3, -5, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 345.88283d0, 0.01917d0, &
   170, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 260.21409d0, 0.01901d0, &
   171, 0, 0, 2, 0, 0, -2, 0, 0, 0, 2, 0, -1, 0, 180.03353d0, 0.01873d0, &
   172, 0, 0, 2, 0, -3, 0, 0, 0, 0, -2, 0, 0, 0, 190.20751d0, 0.01841d0, &
   173, 0, 0, 1, 0, -2, 0, 0, 0, 0, 2, 0, -1, 0, 301.71105d0, 0.01819d0, &
   174, 0, 0, 2, 0, -3, 0, 0, 0, 0, -2, 0, 1, 0, 190.29786d0, 0.01814d0, &
   175, 0, 5, -7, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 169.26762d0, 0.01810d0, &
   176, 0, 0, 3, -7, 0, 0, 0, 0, 0, -1, 0, 2, -1, 144.43917d0, 0.01806d0, &
   177, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 298.04083d0, 0.01770d0, &
   178, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0, 0, -1, 0, 297.82534d0, 0.01712d0, &
   179, 0, 0, 8, -16, 4, 5, 0, 0, 0, 2, 0, -1, 0, 252.68044d0, 0.01709d0, &
   180, 0, 0, 8, -16, 4, 5, 0, 0, 0, -2, 0, 1, 0, 252.68044d0, 0.01709d0, &
   181, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95.38394d0, 0.01699d0, &
   182, 0, 0, 1, 0, 2, -5, 0, 0, 0, 0, 0, 0, 0, 65.64939d0, 0.01696d0, &
   183, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2, 0, -2, 0, 331.68291d0, 0.01682d0, &
   184, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 271.34254d0, 0.01673d0, &
   185, 0, 0, 2, 0, -1, 0, 0, 0, 0, -2, 0, 0, 0, 264.37885d0, 0.01644d0, &
   186, 0, 0, 2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 270.16044d0, 0.01589d0, &
   187, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0.09994d0, 0.01585d0, &
   188, 0, 7, -7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 179.98587d0, 0.01547d0, &
   189, 0, 0, 2, 0, -1, 0, 0, 0, 0, -2, 0, 1, 0, 264.83212d0, 0.01537d0, &
   190, 0, 3, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 2.82174d0, 0.01502d0, &
   191, 0, 3, -5, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 338.66373d0, 0.01498d0, &
   192, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -3, 180.00086d0, 0.01492d0, &
   193, 0, 3, -3, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 2.66594d0, 0.01487d0, &
   194, 0, 15, -12, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 99.70101d0, 0.01467d0, &
   195, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 82.23299d0, 0.01457d0, &
   196, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, -1, 0, 1.96510d0, 0.01455d0, &
   197, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 182.64615d0, 0.01454d0, &
   198, 0, 6, -8, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 347.12979d0, 0.01441d0, &
   199, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, -2, 0, 359.98542d0, 0.01434d0, &
   200, 0, 0, 3, -4, 0, 0, 0, 0, 0, -2, 0, 1, 0, 328.67574d0, 0.01431d0, &
   201, 0, 0, 3, -4, 0, 0, 0, 0, 0, 1, 0, -1, 0, 331.25359d0, 0.01425d0, &
   202, 0, 3, -5, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 346.49874d0, 0.01424d0, &
   203, 0, 3, -5, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 338.06720d0, 0.01423d0, &
   204, 0, 0, 4, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 301.21730d0, 0.01410d0, &
   205, 0, 0, 0, 0, 2, -5, 0, 0, 0, 2, 0, -1, 0, 178.60854d0, 0.01408d0, &
   206, 0, 0, 0, 0, 2, -5, 0, 0, 0, -2, 0, 1, 0, 178.71268d0, 0.01400d0, &
   207, 0, 2, -1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 257.09235d0, 0.01345d0, &
   208, 0, 0, 2, -2, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0.45156d0, 0.01315d0, &
   209, 0, 0, 2, 0, -4, 0, 0, 0, 0, 2, 0, -1, 0, 19.76044d0, 0.01312d0, &
   210, 0, 0, 2, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 307.24896d0, 0.01299d0, &
   211, 0, 1, -1, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 179.87973d0, 0.01292d0, &
   212, 0, 7, -7, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 180.05152d0, 0.01290d0, &
   213, 0, 0, 2, -4, 0, 0, 0, 0, 0, 2, 0, -1, 0, 303.41117d0, 0.01284d0, &
   214, 0, 23, -25, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 168.18245d0, 0.01249d0, &
   215, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 1, -1, 0.00102d0, 0.01234d0, &
   216, 0, 0, 2, 0, -2, 0, 0, 0, 0, 4, 0, -2, 0, 180.12877d0, 0.01229d0, &
   217, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 180.02010d0, 0.01218d0, &
   218, 0, 0, 1, 0, 0, -1, 0, 0, 0, -2, 0, 1, 0, 0.68492d0, 0.01211d0, &
   219, 0, 2, -1, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 257.14013d0, 0.01203d0, &
   220, 0, 0, 2, -2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 180.27603d0, 0.01201d0, &
   221, 0, 3, -4, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 90.34969d0, 0.01197d0, &
   222, 0, 0, 1, 0, -2, 0, 0, 0, 0, 1, 0, 0, -1, 100.53010d0, 0.01183d0, &
   223, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 2, 0, 284.98981d0, 0.01162d0, &
   224, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, -2, 0, 284.98930d0, 0.01162d0, &
   225, 0, 18, -15, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 283.60558d0, 0.01141d0, &
   226, 0, 18, -17, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 129.41497d0, 0.01137d0, &
   227, 0, 0, 3, -4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 331.23293d0, 0.01126d0, &
   228, 0, 0, 2, 0, -3, 1, 0, 0, 0, 2, 0, -2, 0, 181.29764d0, 0.01126d0, &
   229, 0, 0, 3, -8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 29.79747d0, 0.01121d0, &
   230, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0.44680d0, 0.01120d0, &
   231, 0, 0, 4, -7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 271.44643d0, 0.01099d0, &
   232, 0, 0, 3, -4, 0, 0, 0, 0, 0, -2, 0, 0, 0, 328.29186d0, 0.01099d0, &
   233, 0, 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 285.55600d0, 0.01096d0, &
   234, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, -1, -1, 163.74686d0, 0.01093d0, &
   235, 0, 0, 5, -8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 187.28372d0, 0.01091d0, &
   236, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, -1, 163.74984d0, 0.01090d0, &
   237, 0, 0, 2, -3, 0, 0, 0, 0, 0, -2, 0, 1, 0, 330.92395d0, 0.01085d0, &
   238, 0, 8, -8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 179.96606d0, 0.01070d0, &
   239, 0, 3, -5, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 344.00427d0, 0.01070d0, &
   240, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.88394d0, 0.01051d0, &
   241, 0, 0, 2, 0, -2, 0, 0, 0, 0, 4, 0, -1, 0, 180.12629d0, 0.01038d0, &
   242, 0, 0, 2, -4, 0, 0, 0, 0, 0, -2, 0, 1, 0, 302.55753d0, 0.01024d0, &
   243, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 180.00907d0, 0.01019d0, &
   244, 0, 24, -24, 0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 179.43126d0, 0.01017d0 /
   !
   !  64      Latitude  Perturbations     order 0
   !
   DATA (k, (mplan0(I, l, 2), I = 1, 13), phase0(l, 2), amplip0(l, 2), l = 1, 64) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 180.00071d0, 8.04508d0, &
   2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 276.68007d0, 1.51021d0, &
   3, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 26.54287d0, 0.63037d0, &
   4, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 26.54272d0, 0.63014d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0.00075d0, 0.45586d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 180.00069d0, 0.41571d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -2, 0.00086d0, 0.32622d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0.00072d0, 0.29854d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, -1, 0, 180.00073d0, 0.08350d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -1, 0.00000d0, 0.08042d0, &
   11, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 95.13234d0, 0.07755d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 1, 0, 0.00071d0, 0.07332d0, &
   13, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 275.13217d0, 0.07245d0, &
   14, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 26.54302d0, 0.06965d0, &
   15, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 26.54416d0, 0.06747d0, &
   16, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 272.30597d0, 0.06663d0, &
   17, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -1, -1, 180.11961d0, 0.05223d0, &
   18, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -1, 1, 180.12138d0, 0.05069d0, &
   19, 0, 0, 1, 0, -2, 0, 0, 0, 0, 1, 0, 0, 0, 100.50752d0, 0.04964d0, &
   20, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 95.13226d0, 0.04834d0, &
   21, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 180.00070d0, 0.04638d0, &
   22, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0.05680d0, 0.04216d0, &
   23, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -2, 0.00109d0, 0.03977d0, &
   24, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 163.75018d0, 0.03884d0, &
   25, 0, 5, -7, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 169.28198d0, 0.03712d0, &
   26, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 75.50853d0, 0.03480d0, &
   27, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, -1, 179.98734d0, 0.03094d0, &
   28, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, 0, 1, 1.77997d0, 0.03056d0, &
   29, 0, 4, -5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 76.51440d0, 0.02884d0, &
   30, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 1, 179.98167d0, 0.02862d0, &
   31, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 103.33779d0, 0.02736d0, &
   32, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 333.18883d0, 0.02492d0, &
   33, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 333.23168d0, 0.02491d0, &
   34, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 180.00066d0, 0.02403d0, &
   35, 0, 4, -4, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 180.03860d0, 0.02359d0, &
   36, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, 0, -1, 286.68658d0, 0.02348d0, &
   37, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 0, 1, 286.68417d0, 0.02346d0, &
   38, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 0.36009d0, 0.02254d0, &
   39, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, 0, -1, 0.03123d0, 0.02194d0, &
   40, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 1, 26.54240d0, 0.02181d0, &
   41, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, -1, 26.54244d0, 0.02178d0, &
   42, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, 1, 10.08526d0, 0.02012d0, &
   43, 0, 1, -1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 179.98525d0, 0.02008d0, &
   44, 0, 2, -2, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0.00073d0, 0.01993d0, &
   45, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, -1, 10.08059d0, 0.01917d0, &
   46, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 0, 0.00067d0, 0.01745d0, &
   47, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, -1, 0, 0.00075d0, 0.01647d0, &
   48, 0, 0, 1, 0, 0, -2, 0, 0, 0, 1, 0, 0, 0, 113.92221d0, 0.01638d0, &
   49, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 275.13228d0, 0.01504d0, &
   50, 0, 0, 1, 0, -1, 0, 0, 0, 0, 2, 0, 0, -1, 0.45426d0, 0.01409d0, &
   51, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, 0, -1, 26.54132d0, 0.01408d0, &
   52, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0.00072d0, 0.01404d0, &
   53, 0, 0, 1, 0, 2, -5, 0, 0, 0, 1, 0, 0, 0, 239.05935d0, 0.01398d0, &
   54, 0, 2, -3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 256.57252d0, 0.01388d0, &
   55, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -2, 1, 26.54238d0, 0.01365d0, &
   56, 0, 2, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 103.29809d0, 0.01350d0, &
   57, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 154.18809d0, 0.01283d0, &
   58, 0, 0, 2, 0, -2, 0, 0, 0, 0, -2, 0, 0, 1, 180.91996d0, 0.01230d0, &
   59, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, -1, 26.54252d0, 0.01216d0, &
   60, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 95.13224d0, 0.01209d0, &
   61, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, 1, 26.54282d0, 0.01205d0, &
   62, 0, 8, -12, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 153.47650d0, 0.01107d0, &
   63, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 256.31567d0, 0.01018d0, &
   64, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, -2, 0.00102d0, 0.01005d0 /
   !
   ! 115,      Distance,  Perturbations,     order 0
   !
   DATA (k, (mplan0(I, l, 3), I = 1, 13), phase0(l, 3), amplip0(l, 3), l = 1, 115) /&
   1, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -1, 0, 90.11969d0, 1.05870d0, &
   2, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 116.54311d0, 0.72783d0, &
   3, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 296.54574d0, 0.68256d0, &
   4, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 89.98187d0, 0.59827d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -1, 270.00126d0, 0.45648d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 90.00128d0, 0.45276d0, &
   7, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 280.06924d0, 0.41011d0, &
   8, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, 0, 0, 91.79862d0, 0.20497d0, &
   9, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 116.54222d0, 0.20473d0, &
   10, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 296.54299d0, 0.20367d0, &
   11, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 90.36386d0, 0.16644d0, &
   12, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 1, 0, 194.98833d0, 0.15780d0, &
   13, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, -1, 0, 14.98841d0, 0.15780d0, &
   14, 0, 0, 1, 0, -1, 0, 0, 0, 0, -2, 0, 1, 0, 91.74578d0, 0.15751d0, &
   15, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, -1, 0, 89.97863d0, 0.14450d0, &
   16, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 270.00993d0, 0.13811d0, &
   17, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 116.53978d0, 0.13477d0, &
   18, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 296.54238d0, 0.12671d0, &
   19, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 91.22751d0, 0.12666d0, &
   20, 0, 1, -1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 269.98523d0, 0.12362d0, &
   21, 0, 2, -2, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 269.99692d0, 0.12047d0, &
   22, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 90.01606d0, 0.11998d0, &
   23, 0, 2, -2, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 90.31081d0, 0.11617d0, &
   24, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, -1, 0, 197.11421d0, 0.11256d0, &
   25, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 1, 0, 17.11263d0, 0.11251d0, &
   26, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, 0, 0, 196.69224d0, 0.11226d0, &
   27, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 0, 0, 16.68897d0, 0.11216d0, &
   28, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, -1, 270.00092d0, 0.10689d0, &
   29, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 271.06726d0, 0.10504d0, &
   30, 0, 0, 2, 0, -2, 0, 0, 0, 0, -2, 0, 0, 0, 270.93928d0, 0.10503d0, &
   31, 0, 1, -1, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 269.98452d0, 0.10060d0, &
   32, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 90.10540d0, 0.09932d0, &
   33, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, -1, 90.00096d0, 0.09554d0, &
   34, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1, 270.00061d0, 0.08508d0, &
   35, 0, 4, -4, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 89.99224d0, 0.07945d0, &
   36, 0, 0, 2, 0, -3, 0, 0, 0, 0, 0, 0, -1, 0, 280.16516d0, 0.07725d0, &
   37, 0, 6, -8, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 77.22087d0, 0.07054d0, &
   38, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 1, 0, 256.56163d0, 0.06313d0, &
   39, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -1, 0, 76.55968d0, 0.06312d0, &
   40, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 90.70888d0, 0.06209d0, &
   41, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, -2, 0, 270.11778d0, 0.06090d0, &
   42, 0, 0, 2, 0, -2, 0, 0, 0, 0, -2, 0, 1, 0, 270.86062d0, 0.06082d0, &
   43, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 180.28888d0, 0.06007d0, &
   44, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 223.46741d0, 0.05896d0, &
   45, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 43.26801d0, 0.05857d0, &
   46, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1, 5.13453d0, 0.05765d0, &
   47, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 88.57600d0, 0.05723d0, &
   48, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 185.13459d0, 0.05703d0, &
   49, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 358.71194d0, 0.05507d0, &
   50, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, -1, -1, 270.00043d0, 0.05343d0, &
   51, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 1.43434d0, 0.05247d0, &
   52, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 178.55815d0, 0.05164d0, &
   53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 270.06412d0, 0.05035d0, &
   54, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 270.14292d0, 0.04924d0, &
   55, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 90.39176d0, 0.04919d0, &
   56, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 1.58991d0, 0.04883d0, &
   57, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 0, 1, 0, 89.84062d0, 0.04609d0, &
   58, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0.54004d0, 0.04556d0, &
   59, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0, 1, 0, 88.87113d0, 0.04396d0, &
   60, 0, 0, 1, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 270.55045d0, 0.04366d0, &
   61, 0, 0, 2, 0, -2, 0, 0, 0, 0, 2, 0, 0, 0, 90.12570d0, 0.04350d0, &
   62, 0, 0, 1, 0, -2, 0, 0, 0, 0, -2, 0, 0, 0, 36.02580d0, 0.04322d0, &
   63, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, -1, 0, 266.09561d0, 0.04282d0, &
   64, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 271.22625d0, 0.04277d0, &
   65, 0, 2, -2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 269.97026d0, 0.04270d0, &
   66, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 0, 0, 0, 85.30070d0, 0.04106d0, &
   67, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2, 0, -1, 0, 61.73416d0, 0.03974d0, &
   68, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 116.54324d0, 0.03968d0, &
   69, 0, 0, 1, 0, -1, 0, 0, 0, 0, 2, 0, 0, 0, 270.46442d0, 0.03875d0, &
   70, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 1, -1, 90.00028d0, 0.03855d0, &
   71, 0, 5, -5, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 90.06831d0, 0.03760d0, &
   72, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 203.58505d0, 0.03741d0, &
   73, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90.04409d0, 0.03740d0, &
   74, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 89.53531d0, 0.03726d0, &
   75, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 296.54559d0, 0.03722d0, &
   76, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 24.97221d0, 0.03716d0, &
   77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90.00000d0, 0.03644d0, &
   78, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 323.44138d0, 0.03577d0, &
   79, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 143.43350d0, 0.03574d0, &
   80, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 256.15478d0, 0.03539d0, &
   81, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 0.61150d0, 0.03489d0, &
   82, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1.19151d0, 0.03488d0, &
   83, 0, 3, -3, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 269.99573d0, 0.03452d0, &
   84, 0, 0, 1, 0, -2, 0, 0, 0, 0, -2, 0, 1, 0, 36.01309d0, 0.03426d0, &
   85, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, -1, 0, 31.07972d0, 0.03351d0, &
   86, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 1, 0, 90.39429d0, 0.03280d0, &
   87, 0, 8, -13, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 328.08476d0, 0.03125d0, &
   88, 0, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 230.25486d0, 0.03123d0, &
   89, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 270.14383d0, 0.03116d0, &
   90, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 148.04700d0, 0.03103d0, &
   91, 0, 0, 6, -8, 0, 0, 0, 0, 0, 2, 0, -1, 0, 212.34364d0, 0.03048d0, &
   92, 0, 8, -13, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 327.34466d0, 0.03005d0, &
   93, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 147.21064d0, 0.02975d0, &
   94, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 243.03414d0, 0.02873d0, &
   95, 0, 0, 1, -2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 48.21733d0, 0.02862d0, &
   96, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 63.10623d0, 0.02854d0, &
   97, 0, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, -1, 0, 231.97848d0, 0.02841d0, &
   98, 0, 2, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 181.01795d0, 0.02834d0, &
   99, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 90.00001d0, 0.02809d0, &
   100, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 181.14487d0, 0.02778d0, &
   101, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 1, 0, 211.44985d0, 0.02756d0, &
   102, 0, 1, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 90.00927d0, 0.02747d0, &
   103, 0, 2, -3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 179.99887d0, 0.02704d0, &
   104, 0, 4, -5, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 1.37914d0, 0.02671d0, &
   105, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 270.23382d0, 0.02611d0, &
   106, 0, 0, 1, -2, 0, 0, 0, 0, 0, -2, 0, 1, 0, 50.88279d0, 0.02591d0, &
   107, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, -1, 270.00126d0, 0.02488d0, &
   108, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, -1, 90.00128d0, 0.02469d0, &
   109, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, 0, 0, 280.38232d0, 0.02460d0, &
   110, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -2, 0, 100.09665d0, 0.02341d0, &
   111, 0, 0, 2, 0, -3, 0, 0, 0, 0, -2, 0, 0, 0, 280.17755d0, 0.02206d0, &
   112, 0, 6, -6, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 90.04177d0, 0.02196d0, &
   113, 0, 0, 3, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 83.90501d0, 0.02139d0, &
   114, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 181.90502d0, 0.02110d0, &
   115, 0, 0, 2, 0, -1, 0, 0, 0, 0, -2, 0, 0, 0, 354.46019d0, 0.02010d0 /
   !
   ! 154      Longitude Perturbations     order 1
   !
   DATA (k, (mplan1(I, l, 1), I = 1, 13), phase1(l, 1), amplip1(l, 1), l = 1, 154) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.00000d0, 1.67680d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 0, 180.00000d0, 0.51642d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 180.00000d0, 0.41383d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0.00000d0, 0.37115d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0.00000d0, 0.27560d0, &
   6, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 114.56550d0, 0.25425d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 0, 0.00000d0, 0.07118d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0.00000d0, 0.06128d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 180.00000d0, 0.04516d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 180.00000d0, 0.04048d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0.00000d0, 0.03747d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, 0, 180.00000d0, 0.03707d0, &
   13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 0, 180.00000d0, 0.03649d0, &
   14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0.00000d0, 0.02438d0, &
   15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -2, 0, 180.00000d0, 0.02165d0, &
   16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0.00000d0, 0.01923d0, &
   17, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 60.66272d0, 0.01443d0, &
   18, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 114.55388d0, 0.01410d0, &
   19, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 114.79544d0, 0.01326d0, &
   20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0.00000d0, 0.01293d0, &
   21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 0.00000d0, 0.01270d0, &
   22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, 0, 180.00000d0, 0.01097d0, &
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0.00000d0, 0.01072d0, &
   24, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2.27463d0, 0.01056d0, &
   25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0.00000d0, 0.00840d0, &
   26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0.00000d0, 0.00734d0, &
   27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -2, 0, 180.00000d0, 0.00686d0, &
   28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -2, 0, 180.00000d0, 0.00631d0, &
   29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0.00000d0, 0.00585d0, &
   30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, -2, 180.00000d0, 0.00539d0, &
   31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, 0, 0, 180.00000d0, 0.00469d0, &
   32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 1, 0, 180.00000d0, 0.00378d0, &
   33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, -2, 0.00000d0, 0.00362d0, &
   34, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 17.46813d0, 0.00353d0, &
   35, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 93.28221d0, 0.00320d0, &
   36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 180.00000d0, 0.00317d0, &
   37, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 180.00000d0, 0.00300d0, &
   38, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 114.55967d0, 0.00298d0, &
   39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 2, 0, 180.00000d0, 0.00295d0, &
   40, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 0, 0, 75.69874d0, 0.00293d0, &
   41, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 114.55114d0, 0.00293d0, &
   42, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 114.57998d0, 0.00291d0, &
   43, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 114.57740d0, 0.00280d0, &
   44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 180.00000d0, 0.00270d0, &
   45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 180.00000d0, 0.00256d0, &
   46, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 194.81311d0, 0.00247d0, &
   47, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 59.33774d0, 0.00244d0, &
   48, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 57.64503d0, 0.00235d0, &
   49, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, -1, 0, 1.32715d0, 0.00202d0, &
   50, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 1, 0, 1.32745d0, 0.00202d0, &
   51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -3, -1, 0, 180.00000d0, 0.00186d0, &
   52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 0.00000d0, 0.00169d0, &
   53, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, -1, 0, 2.12840d0, 0.00168d0, &
   54, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 1, 0, 2.12715d0, 0.00167d0, &
   55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, -1, 0, 0.00000d0, 0.00160d0, &
   56, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18.08492d0, 0.00157d0, &
   57, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -2, 0, 59.31457d0, 0.00157d0, &
   58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -2, -1, 0, 180.00000d0, 0.00155d0, &
   59, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0.00000d0, 0.00155d0, &
   60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -2, 0, 180.00000d0, 0.00148d0, &
   61, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27.88777d0, 0.00143d0, &
   62, 0, 0, 8, -15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 240.85911d0, 0.00143d0, &
   63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0.00000d0, 0.00141d0, &
   64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0.00000d0, 0.00138d0, &
   65, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 1, 0, 79.70215d0, 0.00131d0, &
   66, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -1, 0, 79.70328d0, 0.00131d0, &
   67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0.00000d0, 0.00128d0, &
   68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -2, 0, 0.00000d0, 0.00126d0, &
   69, 0, 12, -8, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 321.11523d0, 0.00122d0, &
   70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -3, 0, 180.00000d0, 0.00120d0, &
   71, 0, 20, -21, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 18.18419d0, 0.00118d0, &
   72, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 43.20942d0, 0.00108d0, &
   73, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, 0, 0, 1.97140d0, 0.00107d0, &
   74, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 2, 0.00000d0, 0.00107d0, &
   75, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 0, 0, 1.97038d0, 0.00107d0, &
   76, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 180.00000d0, 0.00104d0, &
   77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0.00000d0, 0.00097d0, &
   78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 2, 0.00000d0, 0.00096d0, &
   79, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 114.55394d0, 0.00096d0, &
   80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0.00000d0, 0.00094d0, &
   81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, -2, 0.00000d0, 0.00094d0, &
   82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, -2, 0, 0.00000d0, 0.00090d0, &
   83, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 114.79537d0, 0.00090d0, &
   84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 180.00000d0, 0.00088d0, &
   85, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -2, -2, 0, 180.00000d0, 0.00079d0, &
   86, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 59.60862d0, 0.00079d0, &
   87, 0, 26, -29, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 238.39622d0, 0.00079d0, &
   88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -2, 180.00000d0, 0.00076d0, &
   89, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 180.00000d0, 0.00076d0, &
   90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0.00000d0, 0.00076d0, &
   91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0.00000d0, 0.00073d0, &
   92, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 0, 0.00000d0, 0.00073d0, &
   93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, 1, 0, 180.00000d0, 0.00071d0, &
   94, 0, 20, -20, 0, 0, 0, 0, 0, 0, -1, 0, 1, -1, 297.35671d0, 0.00071d0, &
   95, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, -1, 0, 180.00000d0, 0.00069d0, &
   96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 180.00000d0, 0.00066d0, &
   97, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 60.17999d0, 0.00065d0, &
   98, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 334.24847d0, 0.00063d0, &
   99, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 19.07630d0, 0.00061d0, &
   100, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 43.40532d0, 0.00059d0, &
   101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, -1, 0, 0.00000d0, 0.00058d0, &
   102, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 14.20298d0, 0.00057d0, &
   103, 0, 0, 14, -23, 0, 0, 0, 0, 0, 2, 0, -2, 0, 187.22667d0, 0.00057d0, &
   104, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 43.10853d0, 0.00057d0, &
   105, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 84.42824d0, 0.00055d0, &
   106, 0, 8, -13, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 84.37067d0, 0.00055d0, &
   107, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 194.81311d0, 0.00053d0, &
   108, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1, 194.81311d0, 0.00053d0, &
   109, 0, 0, 6, -8, 0, 0, 0, 0, 0, 2, 0, -1, 0, 44.63031d0, 0.00053d0, &
   110, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 194.81311d0, 0.00052d0, &
   111, 0, 15, -13, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 108.32836d0, 0.00052d0, &
   112, 0, 0, 2, 0, -3, 0, 0, 0, 0, 0, 0, -1, 0, 239.68341d0, 0.00052d0, &
   113, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 10.82813d0, 0.00051d0, &
   114, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, 294.56550d0, 0.00051d0, &
   115, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, 294.56550d0, 0.00051d0, &
   116, 0, 6, -8, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 16.69061d0, 0.00051d0, &
   117, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 19.76925d0, 0.00050d0, &
   118, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 96.01633d0, 0.00049d0, &
   119, 0, 3, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 169.13282d0, 0.00049d0, &
   120, 0, 8, -13, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 96.66542d0, 0.00048d0, &
   121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -3, 0, 180.00000d0, 0.00045d0, &
   122, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 357.99442d0, 0.00045d0, &
   123, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 16.37097d0, 0.00043d0, &
   124, 0, 6, -8, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 359.48088d0, 0.00042d0, &
   125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, -1, 0, 0.00000d0, 0.00041d0, &
   126, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2, 0, -1, 0, 246.46578d0, 0.00041d0, &
   127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0.00000d0, 0.00039d0, &
   128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 180.00000d0, 0.00038d0, &
   129, 0, 3, -4, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 13.81478d0, 0.00038d0, &
   130, 0, 3, -4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 11.16568d0, 0.00038d0, &
   131, 0, 0, 1, 0, 2, -5, 0, 0, 0, 0, 0, 0, 0, 162.46050d0, 0.00037d0, &
   132, 0, 0, 1, 0, -2, 0, 0, 0, 0, -2, 0, 1, 0, 44.09499d0, 0.00036d0, &
   133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0.00000d0, 0.00035d0, &
   134, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, 114.56550d0, 0.00035d0, &
   135, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 114.56550d0, 0.00035d0, &
   136, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 35.03037d0, 0.00035d0, &
   137, 0, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, -1, 0, 57.99304d0, 0.00035d0, &
   138, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 0, 1, 0, 307.21040d0, 0.00035d0, &
   139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0.00000d0, 0.00034d0, &
   140, 0, 0, 6, -8, 0, 0, 0, 0, 0, 2, 0, -2, 0, 35.56634d0, 0.00034d0, &
   141, 0, 2, -3, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 16.76641d0, 0.00034d0, &
   142, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0.00000d0, 0.00033d0, &
   143, 0, 0, 7, -13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 274.01288d0, 0.00033d0, &
   144, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, -2, 180.00000d0, 0.00032d0, &
   145, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 2, 0, 180.00000d0, 0.00032d0, &
   146, 0, 0, 0, 0, 2, -5, 0, 0, 0, 2, 0, -1, 0, 276.56701d0, 0.00032d0, &
   147, 0, 0, 0, 0, 2, -5, 0, 0, 0, -2, 0, 1, 0, 276.71860d0, 0.00032d0, &
   148, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, -1, 0, 307.69915d0, 0.00032d0, &
   149, 0, 2, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 23.49623d0, 0.00032d0, &
   150, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 180.00000d0, 0.00031d0, &
   151, 0, 3, -4, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 0.36084d0, 0.00031d0, &
   152, 0, 0, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62.19489d0, 0.00031d0, &
   153, 0, 8, -13, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 92.98384d0, 0.00031d0, &
   154, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, 0.00000d0, 0.00030d0 /
   !
   !  64      Latitude  Perturbations     order 1
   !
   DATA (k, (mplan1(I, l, 2), I = 1, 13), phase1(l, 2), amplip1(l, 2), l = 1, 64) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, -1, 180.00000d0, 0.07430d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, -1, 0.00000d0, 0.03043d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 1, 180.00000d0, 0.02229d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 1, 180.00000d0, 0.01999d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, -1, 180.00000d0, 0.01869d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0.00000d0, 0.01696d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0.00000d0, 0.01623d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0.00000d0, 0.01419d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0.00000d0, 0.01338d0, &
   10, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 14.44116d0, 0.01304d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0.00000d0, 0.01279d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0.00000d0, 0.01215d0, &
   13, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 114.50884d0, 0.01126d0, &
   14, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 114.62374d0, 0.01123d0, &
   15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, -1, 180.00000d0, 0.00546d0, &
   16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, -1, 180.00000d0, 0.00443d0, &
   17, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.00000d0, 0.00342d0, &
   18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 1, 0.00000d0, 0.00330d0, &
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0.00000d0, 0.00318d0, &
   20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0.00000d0, 0.00295d0, &
   21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 1, 180.00000d0, 0.00285d0, &
   22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, -1, 0.00000d0, 0.00207d0, &
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 180.00000d0, 0.00202d0, &
   24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 180.00000d0, 0.00202d0, &
   25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, -1, 0.00000d0, 0.00200d0, &
   26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, -1, 0.00000d0, 0.00198d0, &
   27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 1, 180.00000d0, 0.00193d0, &
   28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -2, -1, 180.00000d0, 0.00164d0, &
   29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0.00000d0, 0.00161d0, &
   30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, 1, 180.00000d0, 0.00158d0, &
   31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, -1, 180.00000d0, 0.00149d0, &
   32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, -1, 180.00000d0, 0.00135d0, &
   33, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 114.55972d0, 0.00125d0, &
   34, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 114.67748d0, 0.00121d0, &
   35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, 0, -1, 180.00000d0, 0.00104d0, &
   36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, 1, 180.00000d0, 0.00085d0, &
   37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, -1, 0.00000d0, 0.00079d0, &
   38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0.00000d0, 0.00076d0, &
   39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, -1, 0.00000d0, 0.00068d0, &
   40, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 194.81311d0, 0.00066d0, &
   41, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 60.19125d0, 0.00064d0, &
   42, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 60.19125d0, 0.00064d0, &
   43, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 14.81311d0, 0.00062d0, &
   44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0.00000d0, 0.00060d0, &
   45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 1, 0.00000d0, 0.00057d0, &
   46, 0, 5, -6, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 15.60353d0, 0.00057d0, &
   47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 0.00000d0, 0.00057d0, &
   48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -2, 1, 180.00000d0, 0.00054d0, &
   49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, -1, 0.00000d0, 0.00053d0, &
   50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0.00000d0, 0.00053d0, &
   51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 0.00000d0, 0.00048d0, &
   52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, -1, -1, 0.00000d0, 0.00043d0, &
   53, 0, 3, -6, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 56.92022d0, 0.00041d0, &
   54, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 194.81311d0, 0.00041d0, &
   55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, 0, 1, 180.00000d0, 0.00040d0, &
   56, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 1, 114.56550d0, 0.00038d0, &
   57, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, -1, 114.56550d0, 0.00038d0, &
   58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 1, 0.00000d0, 0.00037d0, &
   59, 0, 4, -5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 352.67276d0, 0.00037d0, &
   60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 1, -1, 180.00000d0, 0.00035d0, &
   61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 2, -1, 180.00000d0, 0.00032d0, &
   62, 0, 0, 1, 0, -2, 0, 0, 0, 0, 1, 0, 0, 0, 333.54238d0, 0.00032d0, &
   63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 2, 1, 180.00000d0, 0.00031d0, &
   64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 0.00000d0, 0.00031d0 /
   !
   !  69      Distance  Perturbations     order 1
   !
   DATA (k, (mplan1(I, l, 3), I = 1, 13), phase1(l, 3), amplip1(l, 3), l = 1, 69) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 90.00000d0, 0.51395d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 0, 90.00000d0, 0.38245d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 90.00000d0, 0.32654d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 270.00000d0, 0.26396d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 270.00000d0, 0.12302d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 270.00000d0, 0.07754d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 0, 270.00000d0, 0.06068d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 90.00000d0, 0.04970d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 90.00000d0, 0.04194d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 0, 90.00000d0, 0.03222d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -2, 0, 270.00000d0, 0.02529d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, 0, 90.00000d0, 0.02490d0, &
   13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 90.00000d0, 0.01764d0, &
   14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 270.00000d0, 0.01449d0, &
   15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 271.81691d0, 0.01356d0, &
   16, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 204.55406d0, 0.01302d0, &
   17, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24.79524d0, 0.01225d0, &
   18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 270.00000d0, 0.01186d0, &
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 90.00000d0, 0.01066d0, &
   20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, 0, 90.00000d0, 0.00993d0, &
   21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 270.00000d0, 0.00658d0, &
   22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 270.00000d0, 0.00633d0, &
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 270.00000d0, 0.00587d0, &
   24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 270.00000d0, 0.00536d0, &
   25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -2, 0, 90.00000d0, 0.00476d0, &
   26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, 0, 0, 90.00000d0, 0.00394d0, &
   27, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 24.56550d0, 0.00364d0, &
   28, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 204.56550d0, 0.00364d0, &
   29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 1, 0, 90.00000d0, 0.00331d0, &
   30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 90.00000d0, 0.00310d0, &
   31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 270.00000d0, 0.00246d0, &
   32, 0, 18, -16, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 204.58914d0, 0.00241d0, &
   33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 90.00000d0, 0.00235d0, &
   34, 0, 18, -16, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 24.57709d0, 0.00226d0, &
   35, 0, 0, 2, 0, -3, 0, 0, 0, 0, 2, 0, -1, 0, 329.33777d0, 0.00225d0, &
   36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 270.00000d0, 0.00214d0, &
   37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 2, 0, 90.00000d0, 0.00213d0, &
   38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90.00000d0, 0.00207d0, &
   39, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, 1, 0, 271.29785d0, 0.00178d0, &
   40, 0, 0, 4, -8, 3, 0, 0, 0, 0, 0, 0, -1, 0, 91.29765d0, 0.00178d0, &
   41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -2, 0, 270.00000d0, 0.00173d0, &
   42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, -2, 270.00000d0, 0.00165d0, &
   43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, -1, 0, 270.00000d0, 0.00146d0, &
   44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -2, -1, 0, 90.00000d0, 0.00140d0, &
   45, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 0, 0, 91.97023d0, 0.00134d0, &
   46, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, 0, 0, 271.97125d0, 0.00134d0, &
   47, 0, 0, 4, -8, 3, 0, 0, 0, 0, 2, 0, -1, 0, 272.13364d0, 0.00132d0, &
   48, 0, 0, 4, -8, 3, 0, 0, 0, 0, -2, 0, 1, 0, 92.13235d0, 0.00131d0, &
   49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 270.00000d0, 0.00126d0, &
   50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -3, 0, 270.00000d0, 0.00125d0, &
   51, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, 1, 0, 349.08218d0, 0.00120d0, &
   52, 0, 0, 0, 0, 2, -5, 0, 0, 0, 0, 0, -1, 0, 169.08218d0, 0.00120d0, &
   53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 270.00000d0, 0.00114d0, &
   54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -3, -1, 0, 90.00000d0, 0.00112d0, &
   55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 90.00000d0, 0.00106d0, &
   56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 270.00000d0, 0.00089d0, &
   57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 270.00000d0, 0.00084d0, &
   58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 2, 270.00000d0, 0.00084d0, &
   59, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, -2, 270.00000d0, 0.00081d0, &
   60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 270.00000d0, 0.00074d0, &
   61, 0, 10, -3, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 149.60900d0, 0.00073d0, &
   62, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 204.45391d0, 0.00072d0, &
   63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 90.00000d0, 0.00069d0, &
   64, 0, 18, -16, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24.68404d0, 0.00068d0, &
   65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 90.00000d0, 0.00064d0, &
   66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, -1, 0, 270.00000d0, 0.00064d0, &
   67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 90.00000d0, 0.00062d0, &
   68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 270.00000d0, 0.00061d0, &
   69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, -2, 0, 270.00000d0, 0.00060d0 /
   !
   !  25      Longitude Perturbations     order 2
   !
   DATA (k, (mplan2(I, l, 1), I = 1, 13), phase2(l, 1), amplip2(l, 1), l = 1, 25) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.00000d0, 0.00487d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 0, 180.00000d0, 0.00150d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 180.00000d0, 0.00120d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0.00000d0, 0.00108d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0.00000d0, 0.00080d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 0, 0.00000d0, 0.00021d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0.00000d0, 0.00018d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 180.00000d0, 0.00013d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 180.00000d0, 0.00012d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 0, 180.00000d0, 0.00011d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, 0, 180.00000d0, 0.00011d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0.00000d0, 0.00011d0, &
   13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0.00000d0, 0.00007d0, &
   14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -2, 0, 180.00000d0, 0.00006d0, &
   15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0.00000d0, 0.00006d0, &
   16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 0.00000d0, 0.00004d0, &
   17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0.00000d0, 0.00004d0, &
   18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, 0, 180.00000d0, 0.00003d0, &
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0.00000d0, 0.00003d0, &
   20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -2, 0, 180.00000d0, 0.00002d0, &
   21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -2, 0, 180.00000d0, 0.00002d0, &
   22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, -2, 180.00000d0, 0.00002d0, &
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0.00000d0, 0.00002d0, &
   24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0.00000d0, 0.00002d0, &
   25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0.00000d0, 0.00002d0 /
   !
   !  12      Latitude  Perturbations     order 2
   !
   DATA (k, (mplan2(I, l, 2), I = 1, 13), phase2(l, 2), amplip2(l, 2), l = 1, 12) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, -1, 180.00000d0, 0.00022d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, -1, 0.00000d0, 0.00009d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 1, 180.00000d0, 0.00006d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 1, 180.00000d0, 0.00006d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, -1, 180.00000d0, 0.00005d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0.00000d0, 0.00005d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0.00000d0, 0.00005d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0.00000d0, 0.00004d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0.00000d0, 0.00004d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0.00000d0, 0.00004d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0.00000d0, 0.00004d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, -1, 180.00000d0, 0.00002d0 /
   !
   !  19      Distance  Perturbations     order 2
   !
   DATA (k, (mplan2(I, l, 3), I = 1, 13), phase2(l, 3), amplip2(l, 3), l = 1, 19) /&
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 90.00000d0, 0.00149d0, &
   2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -1, 0, 90.00000d0, 0.00111d0, &
   3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 90.00000d0, 0.00095d0, &
   4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 270.00000d0, 0.00077d0, &
   5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 270.00000d0, 0.00036d0, &
   6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 270.00000d0, 0.00023d0, &
   7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, 0, 270.00000d0, 0.00018d0, &
   8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 90.00000d0, 0.00014d0, &
   9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 90.00000d0, 0.00012d0, &
   10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1, 0, 90.00000d0, 0.00009d0, &
   11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, -2, 0, 270.00000d0, 0.00007d0, &
   12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, 0, 90.00000d0, 0.00007d0, &
   13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 90.00000d0, 0.00005d0, &
   14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 270.00000d0, 0.00004d0, &
   15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 270.00000d0, 0.00004d0, &
   16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, 0, 90.00000d0, 0.00003d0, &
   17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 270.00000d0, 0.00003d0, &
   18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 90.00000d0, 0.00003d0, &
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 270.00000d0, 0.00002d0 /



   !
   ! computation of the light travel time for the aberration correction
   !
   deltat = 0d0
   if(iopt == 2) then ! apparent position computed
      deltat = dist_moon(xjd)/vlight/days ! light time in days
   endif
   tt = ((xjd-xjd2000) - deltat)/century ! time argument in julian century

   tt2 = tt * tt
   tt3 = tt2 * tt
   tt4 = tt2 * tt2


   !
   ! Arguments at xjd for the main problem and the planetary perturbations
   !

   do id = 1, 4 ! loop on Delaunay arguments
      darg(id) = delau0(id) + (delau1(id) * tt + delau2(id) * tt2 + delau3(id) * tt3 + delau4(id) * tt4)/3600d0
   enddo

   do ip = 1, 8 ! loop on the eight planets
      parg(ip) = planet0(ip) + planet1(ip) * tt/3600d0
   enddo

   !
   !  mean longitude of the Moon equinox of date (gamma ICRS)
   !
   xlong = (xlmoon(0) + (xlmoon(1) + cpreces(1)) * tt + (xlmoon(2) + cpreces(2)) * tt2 + (xlmoon(3) + cpreces(3)) * tt3 + (xlmoon(4) + cpreces(4)) * tt4)/3600d0


   ! 13 arguments for planetary perturbations
   !

   parg(9) = xlong
   parg(10:13) = darg(1:4)

   darg = modulo(darg, 360d0)
   parg = modulo(parg, 360d0)

   !
   !  Perturbations in longitude, latitude, distance
   !

   !
   ! Loop truncations according to accuracy level
   ! iprec = 2  :: full precision
   ! iprec = 1  ::  ~ 0"5 and 1 km
   ! iprec = 0  !!  ~ 5"  and 10 km
   ! Amplitudes selected with a margin factor of 5
   !

   x = 0d0
   do k = 1, 3 ! longitude, latitude, distance (deg, deg, km)

      lmax = idmain(k, iprec) ! determined by the accuracy level

      do l = 1, lmax ! main problem (Kepler + Sun)
         arg = 0d0
         do i = 1, 4
            arg = arg + mcent(i, l, k) * darg(i)
         enddo
         arg = arg * degrad
         if (k == 3) then
            x(k) = x(k) + ampli(l, k) * cos(arg)
         else
            x(k) = x(k) + ampli(l, k) * sin(arg)
         endif
      enddo

      lmax = idpla0(k, iprec) ! determined by the accuracy level

      do l = 1, lmax ! planetary pert. order 0
         arg = phase0(l, k)
         do i = 1, 13
            arg = arg + mplan0(i, l, k) * parg(i)
         enddo
         x(k) = x(k) + amplip0(l, k) * sin(arg * degrad)
      enddo

      lmax = idpla1(k, iprec) ! determined by the accuracy level

      y = 0d0
      do l = 1, lmax ! planetary pert. order 1
         arg = phase1(l, k)
         do i = 1, 13
            arg = arg + mplan1(i, l, k) * parg(i)
         enddo
         y(k) = y(k) + amplip1(l, k) * sin(arg * degrad)
      enddo
      x = x + tt * y


      lmax = idpla2(k, iprec) ! determined by the accuracy level

      y = 0d0
      do l = 1, lmax ! planetary pert. order 2
         arg = phase2(l, k)
         do i = 1, 13
            arg = arg + mplan2(i, l, k) * parg(i)
         enddo
         y(k) = y(k) + amplip2(l, k) * sin(arg * degrad)
      enddo
      x = x + tt2 * y
   enddo ! end of loop on coordinates


   xlong = modulo(xlong + x(1)/3600d0, round) ! mean long + elliptic + perturbations
   xlat = x(2)/3600d0
   rho = x(3)

   !
   ! Apparent position and equatorial coordinates
   !

   call dpsideps(xjd, dpsi, deps) ! nutation in longitude and obliquity
   if(iopt == 1) then
      epsilon = obliquity(xjd)
   else
      epsilon = obliquity(xjd) + deps/3600d0
      xlong = xlong + dpsi/3600d0 ! nutation included
   endif


   call sphecar(rho, xlong, xlat, posm) ! position vector of the moon in ecliptic frame

   !
   ! Equatorial rectangular coordinates (ICRS origin in the equator)
   !

   call rota(posm, -epsilon, 1, pos)

   call rota(pos, phid, 3, pos) ! this must be improved with NRO if t /= J2000

   !
   ! Equatorial spherical coordinates
   !
   call carsphe(pos, rho, ra, dec)

   return
end subroutine Moon
function dist_moon(xjd)
   !*****************************************************************
   ! Earth Moon distance at julian day xjd
   ! dist_moon in km
   ! Low precision (~ 100km) used to determine the light time
   ! to compute the aberration
   !
   ! F. Mignard  December 2005 & November 2016
   !
   ! INPUT
   !   xjd   : julian day
   !
   ! OUTPUT
   ! dist_moon : Earth-Moon distance in km
   !******************************************************************
   implicit none

   real(kind = dp), intent(in) :: xjd
   real(kind = dp) :: dist_moon

   real(kind = dp) :: tt, rho
   integer :: k
   integer, parameter          :: nterm = 15
   real(kind = dp), parameter :: rhom = 385000.53d0
!   real(kind = dp), dimension(4), parameter :: ampli = (/-20905d0, -3699d0, -2956d0, -569.9d0/) ! km
!   real(kind = dp), dimension(4), parameter :: phase = (/134.9d0, 100.7d0, 235.7d0, 269.9d0/) ! deg
!   real(kind = dp), dimension(4), parameter :: freq = (/477198.868d0, 413335.355d0, 890534.223d0, 954397.735d0/)! deg/cy

   real(kind = dp), dimension(nterm), parameter :: ampli =&
   (/ -20905.3d0, -3699.1d0, -2955.9d0, -569.9d0, 246.1d0, -204.6d0, -170.7d0, -152.1d0, -129.6d0, 108.7d0,104.8d0, 79.7d0, 48.9d0, -34.8d0, 30.8d0/) ! km
   real(kind = dp), dimension(nterm), parameter :: phase =&
   (/ 134.96d0, 100.74d0, 235.70d0, 269.93d0, 325.77d0, 238.17d0, 10.66d0, 103.21d0, 222.57d0, 297.85d0, 132.49d0, 308.42d0, 357.53d0, 336.44d0, 233.23d0/) ! deg    
   real(kind = dp), dimension(nterm), parameter :: freq =&
  (/ 477198.868d0, 413335.355d0,  890534.223d0, 954397.735d0,  -63863.512d0,  854535.173d0,&
    1367733.090d0, 377336.305d0, -441199.817d0, 445267.112d0,  513197.918d0, -489205.167d0, 35999.050d0, 1303869.578d0, 926533.273d0/)! deg/cy
   
   
   tt = (xjd-xjd2000)/century
   rho = rhom
   do k = 1, nterm
      rho = rho + ampli(k) * cos((freq(k) * tt + phase(k)) * degrad)
   enddo
   dist_moon = rho
   return
end function dist_moon
subroutine vit_moon(datjd, iframe, vit)
   !
   !     Velocity of the Moon  with the origin at the
   !     Earth barycenter.
   !     Mean equator (or mean ecliptic) and mean equinox of date
   !
   !
   !     F. Mignard March 2003
   !
   !***************************************************************
   ! INPUT
   !     datjd   : date in julian days
   !     iframe  : 0   =  Equatorial frame, 1 = ecliptic frame
   !
   ! OUTPUT
   !     vit     : array(3) with the velocity vector in km/s
   !***************************************************************

   !
   !
   implicit none

   real(kind = dp) datjd
   real(kind = dp) vit(3)
   integer iframe
   real(kind = dp) pos1(3), pos2(3), pos3(3), pos4(3), pos(3)
   real(kind = dp) xlong, xlat, rho, ra, dec


   integer, parameter :: iprec = 1 ! low precision orbit for the Moon
   integer, parameter :: iopt = 1 ! geometric Moon


   real(kind = dp), parameter :: dt = 0.05d0 ! step in days for derivation optimized for period of ~ 1 month


   call moon(datjd-1d0 * dt, iprec, iopt, xlong, xlat, rho, ra, dec, pos1)
   call moon(datjd + 1d0 * dt, iprec, iopt, xlong, xlat, rho, ra, dec, pos2)
   call moon(datjd-2d0 * dt, iprec, iopt, xlong, xlat, rho, ra, dec, pos3)
   call moon(datjd + 2d0 * dt, iprec, iopt, xlong, xlat, rho, ra, dec, pos4)


   vit = ((pos3 - pos4) + 8 * (pos2 - pos1))/(12 * dt)/days !speed in km/s in equatorial frame

   if (iframe == 1) then
      call equtoecl(datjd, vit, vit) ! passage to ecliptic frame
   endif

   return
end subroutine vit_moon

!********************************************************
subroutine baryc_l(xjd, xyz, xyzp)
   !********************************************************
   !
   !     Vector Barycenter to  Sun at the julian date xjd
   !     xyz(i)  in au,  xyzp in au/day
   !     ecliptic frame at J2000.
   !     Low precision
   !
   !     Source ! from VSOP Series provided by P. Bretagnon June 2000
   !
   !     F. Mignard OCA/Cassiopee
   !
   !    *INPUT
   !     xjd     : julian day
   !
   !    *OUTPUT
   !     xyz     : Sun from SS barycenter in au in ecliptic frame J2000
   !     xyzp    : Velocity of the Sun in au/day
   !
   !     Precision : 1.0e-5 au      = 1500 km
   !                 2.0e-7 au/day  = 0.3  m/s
   !********************************************************
   !
   implicit none

   real(kind = dp), intent(in) :: xjd
   real(kind = dp), dimension(3), intent(out) :: xyz, xyzp
   real(kind = dp) :: tab(3, 30, 3)
   integer :: kf(3) = (/30, 30, 7/)
   real(kind = dp) :: t, aa, bb, arg, amp
   integer :: j, k
   real(kind = dp), parameter :: millenia = 365250d0



   data tab/&
   0.00495672739d0, 3.74107356792d0, 529.69096509460d0, &
   0.00271802376d0, 4.01601149775d0, 213.29909543800d0, &
   0.00155435675d0, 2.17052050061d0, 38.13303563780d0, &
   0.00083792997d0, 2.33967985523d0, 74.78159856730d0, &
   0.00029374249d0, 0.00000000000d0, 0.00000000000d0, &
   0.00012013079d0, 4.09073224903d0, 1059.38193018920d0, &
   0.00007577257d0, 3.24151897354d0, 426.59819087600d0, &
   0.00001941380d0, 1.01219474101d0, 206.18554843720d0, &
   0.00001940649d0, 4.79769963661d0, 149.56319713460d0, &
   0.00001888831d0, 3.89252804366d0, 220.41264243880d0, &
   0.00001434208d0, 3.86895363775d0, 522.57741809380d0, &
   0.00001406367d0, 0.47598335150d0, 536.80451209540d0, &
   0.00001185835d0, 0.77770585045d0, 76.26607127560d0, &
   0.00000813685d0, 3.25483611884d0, 36.64856292950d0, &
   0.00000767074d0, 4.22743731914d0, 39.61750834610d0, &
   0.00000624814d0, 0.27936466811d0, 73.29712585900d0, &
   0.00000436640d0, 4.44044655092d0, 1589.07289528380d0, &
   0.00000379145d0, 5.15640874752d0, 7.11354700080d0, &
   0.00000315393d0, 6.15699854629d0, 419.48464387520d0, &
   0.00000308784d0, 2.49456658747d0, 639.89728631400d0, &
   0.00000278795d0, 4.93431870348d0, 110.20632121940d0, &
   0.00000303993d0, 4.89507841707d0, 6283.07584999140d0, &
   0.00000227188d0, 5.27839813806d0, 103.09277421860d0, &
   0.00000216162d0, 5.80298032120d0, 316.39186965660d0, &
   0.00000176764d0, 0.03416525046d0, 10213.28554621100d0, &
   0.00000135792d0, 2.00151020964d0, 1.48447270830d0, &
   0.00000116993d0, 2.42475255571d0, 632.78373931320d0, &
   0.00000105413d0, 3.12332213850d0, 433.71173787680d0, &
   0.00000097988d0, 3.02626461372d0, 1052.26838318840d0, &
   0.00000109101d0, 3.15781282608d0, 1162.47470440780d0, &

   0.00495536218d0, 2.17046712634d0, 529.69096509460d0, &
   0.00272185821d0, 2.44443364925d0, 213.29909543800d0, &
   0.00155444313d0, 0.59927010840d0, 38.13303563780d0, &
   0.00083755792d0, 0.76880164710d0, 74.78159856730d0, &
   0.00033869535d0, 0.00000000000d0, 0.00000000000d0, &
   0.00012011827d0, 2.52003147880d0, 1059.38193018920d0, &
   0.00007585830d0, 1.66995483217d0, 426.59819087600d0, &
   0.00001963743d0, 5.70773655842d0, 206.18554843720d0, &
   0.00001891503d0, 2.32096821003d0, 220.41264243880d0, &
   0.00001940124d0, 3.22686130461d0, 149.56319713460d0, &
   0.00001436841d0, 2.30161968078d0, 522.57741809380d0, &
   0.00001405975d0, 5.18858607879d0, 536.80451209540d0, &
   0.00001185515d0, 5.48969329104d0, 76.26607127560d0, &
   0.00000813077d0, 1.68393442622d0, 36.64856292950d0, &
   0.00000767125d0, 2.65620459324d0, 39.61750834610d0, &
   0.00000628788d0, 4.99295631526d0, 73.29712585900d0, &
   0.00000436632d0, 2.86969820654d0, 1589.07289528380d0, &
   0.00000382844d0, 3.57213982765d0, 7.11354700080d0, &
   0.00000317511d0, 4.53536380695d0, 419.48464387520d0, &
   0.00000309191d0, 0.92301535903d0, 639.89728631400d0, &
   0.00000287366d0, 3.36314089821d0, 110.20632121940d0, &
   0.00000304013d0, 3.32425157103d0, 6283.07584999140d0, &
   0.00000269924d0, 0.29178785093d0, 103.09277421860d0, &
   0.00000213445d0, 4.22099738237d0, 316.39186965660d0, &
   0.00000177041d0, 4.74733135300d0, 10213.28554621100d0, &
   0.00000138577d0, 0.43043981485d0, 1.48447270830d0, &
   0.00000112761d0, 0.85382170184d0, 632.78373931320d0, &
   0.00000105538d0, 1.55181188435d0, 433.71173787680d0, &
   0.00000098007d0, 1.45965911177d0, 1052.26838318840d0, &
   0.00000109014d0, 1.58735183284d0, 1162.47470440780d0, &

   0.00011810648d0, 0.46078690233d0, 213.29909543800d0, &
   0.00011277700d0, 0.41689943638d0, 529.69096509460d0, &
   0.00004802048d0, 4.58264723370d0, 38.13303563780d0, &
   0.00001131046d0, 5.75877139035d0, 74.78159856730d0, &
   0.00001152656d0, 3.14159265359d0, 0.00000000000d0, &
   0.00000329820d0, 5.97879747107d0, 426.59819087600d0, &
   0.00000273335d0, 0.76652182727d0, 1059.38193018920d0, &
   69 * 0.d0/

   t = (xjd-xjd2000)/millenia
   do j = 1, 3
      aa = 0d0
      bb = 0d0
      do k = 1, kf(j)
         arg = tab(2, k, j) + tab(3, k, j) * t
         amp = tab(1, k, j)
         aa = aa + amp * cos(arg)
         bb = bb - tab(3, k, j) * amp * sin(arg)
      enddo
      xyz(j) = aa ! position in au
      xyzp(j) = bb/millenia ! velocity in au/day
   enddo
   return
end subroutine baryc_l
!********************************************************
subroutine baryc_m(xjd, xyz, xyzp)
   !********************************************************
   !
   !     Vector Barycenter to  Sun at the julian date xjd
   !     xyz(i)  in au,  xyzp in au/day
   !     ecliptic frame at J2000.
   !     Medium precision
   !
   !     Source ! from VSOP Series provided by P. Bretagnon June 2000
   !
   !     F. Mignard OCA/Cassiopee
   !
   !    *INPUT
   !     xjd     : julian day
   !
   !    *OUTPUT
   !     xyz     : Sun from SS barycenter in au in ecliptic frame J2000
   !     xyzp    : Velocity of the Sun in au/day
   !
   !     Precision : 6.0e-7 au      ~ 100 km
   !                 1.0e-8 au/day  ~ 20 mm/s
   !
   !********************************************************
   !
   implicit none

   integer, parameter :: ns0m = 151, ns1m = 40, ns2m = 9
   integer, parameter :: ns0x = 151, ns1x = 40, ns2x = 9
   integer, parameter :: ns0y = 151, ns1y = 40, ns2y = 9
   integer, parameter :: ns0z = 51, ns1z = 14, ns2z = 2

   real(kind = dp), parameter :: poseps = 1.0d-7 ! 1/5 accuracy in position in au
   real(kind = dp), parameter :: viteps = 5.0d-7 ! 1/5 accuracy in velocity in au/yr

   real(kind = dp), intent(in) :: xjd
   real(kind = dp), dimension(3), intent(out) :: xyz, xyzp


   !max number of terms for the required accuracy. Determined by trial and errors. Adapted to the above precision.

   integer :: ns0(3) = (/80, 80, 20/)
   integer :: ns1(3) = (/10, 10, 0/)
   integer :: ns2(3) = (/0, 0, 0/)

   !      integer                             :: ns0(3) =(/ns0x,ns0y,ns0z/)
   !      integer                             :: ns1(3) =(/ns1x,ns1y,ns1z/)
   !      integer                             :: ns2(3) =(/0,0,0/)

   real(kind = dp) :: s0(3, ns0m, 3), s1(3, ns1m, 3), s2(3, ns2m, 3)

   real(kind = dp) :: t, t2, aa, bb, arg, amp, cc, freq, ampv
   integer :: i, j, k


   !SSB-to-Sun, T^0, X
   DATA ((S0(I, J, 1), I = 1, 3), J = 1, 10) /&
   0.4956757536410d-02, 0.3741073751789d+01, 0.5296909721118d+00, &
   0.2718490072522d-02, 0.4016011511425d+01, 0.2132990797783d+00, &
   0.1546493974344d-02, 0.2170528330642d+01, 0.3813291813120d-01, &
   0.8366855276341d-03, 0.2339614075294d+01, 0.7478166569050d-01, &
   0.2936777942117d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.1201317439469d-03, 0.4090736353305d+01, 0.1059381944224d+01, &
   0.7578550887230d-04, 0.3241518088140d+01, 0.4265981595566d+00, &
   0.1941787367773d-04, 0.1012202064330d+01, 0.2061856251104d+00, &
   0.1889227765991d-04, 0.3892520416440d+01, 0.2204125344462d+00, &
   0.1937896968613d-04, 0.4797779441161d+01, 0.1495633313810d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 11, 20) /&
   0.1434506110873d-04, 0.3868960697933d+01, 0.5225775174439d+00, &
   0.1406659911580d-04, 0.4759766557397d+00, 0.5368044267797d+00, &
   0.1179022300202d-04, 0.7774961520598d+00, 0.7626583626240d-01, &
   0.8085864460959d-05, 0.3254654471465d+01, 0.3664874755930d-01, &
   0.7622752967615d-05, 0.4227633103489d+01, 0.3961708870310d-01, &
   0.6209171139066d-05, 0.2791828325711d+00, 0.7329749511860d-01, &
   0.4366435633970d-05, 0.4440454875925d+01, 0.1589072916335d+01, &
   0.3792124889348d-05, 0.5156393842356d+01, 0.7113454667900d-02, &
   0.3154548963402d-05, 0.6157005730093d+01, 0.4194847048887d+00, &
   0.3088359882942d-05, 0.2494567553163d+01, 0.6398972393349d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 21, 30) /&
   0.2788440902136d-05, 0.4934318747989d+01, 0.1102062672231d+00, &
   0.3039928456376d-05, 0.4895077702640d+01, 0.6283075850446d+01, &
   0.2272258457679d-05, 0.5278394064764d+01, 0.1030928125552d+00, &
   0.2162007057957d-05, 0.5802978019099d+01, 0.3163918923335d+00, &
   0.1767632855737d-05, 0.3415346595193d-01, 0.1021328554739d+02, &
   0.1349413459362d-05, 0.2001643230755d+01, 0.1484170571900d-02, &
   0.1170141900476d-05, 0.2424750491620d+01, 0.6327837846670d+00, &
   0.1054355266820d-05, 0.3123311487576d+01, 0.4337116142245d+00, &
   0.9800822461610d-06, 0.3026258088130d+01, 0.1052268489556d+01, &
   0.1091203749931d-05, 0.3157811670347d+01, 0.1162474756779d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 31, 40) /&
   0.6960236715913d-06, 0.8219570542313d+00, 0.1066495398892d+01, &
   0.5689257296909d-06, 0.1323052375236d+01, 0.9491756770005d+00, &
   0.6613172135802d-06, 0.2765348881598d+00, 0.8460828644453d+00, &
   0.6277702517571d-06, 0.5794064466382d+01, 0.1480791608091d+00, &
   0.6304884066699d-06, 0.7323555380787d+00, 0.2243449970715d+00, &
   0.4897850467382d-06, 0.3062464235399d+01, 0.3340612434717d+01, &
   0.3759148598786d-06, 0.4588290469664d+01, 0.3516457698740d-01, &
   0.3110520548195d-06, 0.1374299536572d+01, 0.6373574839730d-01, &
   0.3064708359780d-06, 0.4222267485047d+01, 0.1104591729320d-01, &
   0.2856347168241d-06, 0.3714202944973d+01, 0.1510475019529d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 41, 50) /&
   0.2840945514288d-06, 0.2847972875882d+01, 0.4110125927500d-01, &
   0.2378951599405d-06, 0.3762072563388d+01, 0.2275259891141d+00, &
   0.2714229481417d-06, 0.1036049980031d+01, 0.2535050500000d-01, &
   0.2323551717307d-06, 0.4682388599076d+00, 0.8582758298370d-01, &
   0.1881790512219d-06, 0.4790565425418d+01, 0.2118763888447d+01, &
   0.2261353968371d-06, 0.1669144912212d+01, 0.7181332454670d-01, &
   0.2214546389848d-06, 0.3937717281614d+01, 0.2968341143800d-02, &
   0.2184915594933d-06, 0.1129169845099d+00, 0.7775000683430d-01, &
   0.2000164937936d-06, 0.4030009638488d+01, 0.2093666171530d+00, &
   0.1966105136719d-06, 0.8745955786834d+00, 0.2172315424036d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 51, 60) /&
   0.1904742332624d-06, 0.5919743598964d+01, 0.2022531624851d+00, &
   0.1657399705031d-06, 0.2549141484884d+01, 0.7358765972222d+00, &
   0.1574070533987d-06, 0.5277533020230d+01, 0.7429900518901d+00, &
   0.1832261651039d-06, 0.3064688127777d+01, 0.3235053470014d+00, &
   0.1733615346569d-06, 0.3011432799094d+01, 0.1385174140878d+00, &
   0.1549124014496d-06, 0.4005569132359d+01, 0.5154640627760d+00, &
   0.1637044713838d-06, 0.1831375966632d+01, 0.8531963191132d+00, &
   0.1123420082383d-06, 0.1180270407578d+01, 0.1990721704425d+00, &
   0.1083754165740d-06, 0.3414101320863d+00, 0.5439178814476d+00, &
   0.1156638012655d-06, 0.6130479452594d+00, 0.5257585094865d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 61, 70) /&
   0.1142548785134d-06, 0.3724761948846d+01, 0.5336234347371d+00, &
   0.7921463895965d-07, 0.2435425589361d+01, 0.1478866649112d+01, &
   0.7428600285231d-07, 0.3542144398753d+01, 0.2164800718209d+00, &
   0.8323211246747d-07, 0.3525058072354d+01, 0.1692165728891d+01, &
   0.7257595116312d-07, 0.1364299431982d+01, 0.2101180877357d+00, &
   0.7111185833236d-07, 0.2460478875808d+01, 0.4155522422634d+00, &
   0.6868090383716d-07, 0.4397327670704d+01, 0.1173197218910d+00, &
   0.7226419974175d-07, 0.4042647308905d+01, 0.1265567569334d+01, &
   0.6955642383177d-07, 0.2865047906085d+01, 0.9562891316684d+00, &
   0.7492139296331d-07, 0.5014278994215d+01, 0.1422690933580d-01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 71, 80) /&
   0.6598363128857d-07, 0.2376730020492d+01, 0.6470106940028d+00, &
   0.7381147293385d-07, 0.3272990384244d+01, 0.1581959461667d+01, &
   0.6402909624032d-07, 0.5302290955138d+01, 0.9597935788730d-01, &
   0.6237454263857d-07, 0.5444144425332d+01, 0.7084920306520d-01, &
   0.5241198544016d-07, 0.4215359579205d+01, 0.5265099800692d+00, &
   0.5144463853918d-07, 0.1218916689916d+00, 0.5328719641544d+00, &
   0.5868164772299d-07, 0.2369402002213d+01, 0.7871412831580d-01, &
   0.6233195669151d-07, 0.1254922242403d+01, 0.2608790314060d+02, &
   0.6068463791422d-07, 0.5679713760431d+01, 0.1114304132498d+00, &
   0.4359361135065d-07, 0.6097219641646d+00, 0.1375773836557d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 81, 90) /&
   0.4686510366826d-07, 0.4786231041431d+01, 0.1143987543936d+00, &
   0.3758977287225d-07, 0.1167368068139d+01, 0.1596186371003d+01, &
   0.4282051974778d-07, 0.1519471064319d+01, 0.2770348281756d+00, &
   0.5153765386113d-07, 0.1860532322984d+01, 0.2228608264996d+00, &
   0.4575129387188d-07, 0.7632857887158d+00, 0.1465949902372d+00, &
   0.3326844933286d-07, 0.1298219485285d+01, 0.5070101000000d-01, &
   0.3748617450984d-07, 0.1046510321062d+01, 0.4903339079539d+00, &
   0.2816756661499d-07, 0.3434522346190d+01, 0.2991266627620d+00, &
   0.3412750405039d-07, 0.2523766270318d+01, 0.3518164938661d+00, &
   0.2655796761776d-07, 0.2904422260194d+01, 0.6256703299991d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 91, 100) /&
   0.2963597929458d-07, 0.5923900431149d+00, 0.1099462426779d+00, &
   0.2539523734781d-07, 0.4851947722567d+01, 0.1256615170089d+02, &
   0.2283087914139d-07, 0.3400498595496d+01, 0.6681224869435d+01, &
   0.2321309799331d-07, 0.5789099148673d+01, 0.3368040641550d-01, &
   0.2549657649750d-07, 0.3991856479792d-01, 0.1169588211447d+01, &
   0.2290462303977d-07, 0.2788567577052d+01, 0.1045155034888d+01, &
   0.1945398522914d-07, 0.3290896998176d+01, 0.1155361302111d+01, &
   0.1849171512638d-07, 0.2698060129367d+01, 0.4452511715700d-02, &
   0.1647199834254d-07, 0.3016735644085d+01, 0.4408250688924d+00, &
   0.1529530765273d-07, 0.5573043116178d+01, 0.6521991896920d-01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 101, 110) /&
   0.1433199339978d-07, 0.1481192356147d+01, 0.9420622223326d+00, &
   0.1729134193602d-07, 0.1422817538933d+01, 0.2108507877249d+00, &
   0.1716463931346d-07, 0.3469468901855d+01, 0.2157473718317d+00, &
   0.1391206061378d-07, 0.6122436220547d+01, 0.4123712502208d+00, &
   0.1404746661924d-07, 0.1647765641936d+01, 0.4258542984690d-01, &
   0.1410452399455d-07, 0.5989729161964d+01, 0.2258291676434d+00, &
   0.1089828772168d-07, 0.2833705509371d+01, 0.4226656969313d+00, &
   0.1047374564948d-07, 0.5090690007331d+00, 0.3092784376656d+00, &
   0.1358279126532d-07, 0.5128990262836d+01, 0.7923417740620d-01, &
   0.1020456476148d-07, 0.9632772880808d+00, 0.1456308687557d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 111, 120) /&
   0.1033428735328d-07, 0.3223779318418d+01, 0.1795258541446d+01, &
   0.1412435841540d-07, 0.2410271572721d+01, 0.1525316725248d+00, &
   0.9722759371574d-08, 0.2333531395690d+01, 0.8434341241180d-01, &
   0.9657334084704d-08, 0.6199270974168d+01, 0.1272681024002d+01, &
   0.1083641148690d-07, 0.2864222292929d+01, 0.7032915397480d-01, &
   0.1067318403838d-07, 0.5833458866568d+00, 0.2123349582968d+00, &
   0.1062366201976d-07, 0.4307753989494d+01, 0.2142632012598d+00, &
   0.1236364149266d-07, 0.2873917870593d+01, 0.1847279083684d+00, &
   0.1092759489593d-07, 0.2959887266733d+01, 0.1370332435159d+00, &
   0.8912069362899d-08, 0.5141213702562d+01, 0.2648454860559d+01 /

   DATA ((S0(I, J, 1), I = 1, 3), J = 121, 130) /&
   0.9656467707970d-08, 0.4532182462323d+01, 0.4376440768498d+00, &
   0.8098386150135d-08, 0.2268906338379d+01, 0.2880807454688d+00, &
   0.7857714675000d-08, 0.4055544260745d+01, 0.2037373330570d+00, &
   0.7288455940646d-08, 0.5357901655142d+01, 0.1129145838217d+00, &
   0.9450595950552d-08, 0.4264926963939d+01, 0.5272426800584d+00, &
   0.9381718247537d-08, 0.7489366976576d-01, 0.5321392641652d+00, &
   0.7079052646038d-08, 0.1923311052874d+01, 0.6288513220417d+00, &
   0.9259004415344d-08, 0.2970256853438d+01, 0.1606092486742d+00, &
   0.8259801499742d-08, 0.3327056314697d+01, 0.8389694097774d+00, &
   0.6476334355779d-08, 0.2954925505727d+01, 0.2008557621224d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 131, 140) /&
   0.5984021492007d-08, 0.9138753105829d+00, 0.2042657109477d+02, &
   0.5989546863181d-08, 0.3244464082031d+01, 0.2111650433779d+01, &
   0.6233108606023d-08, 0.4995232638403d+00, 0.4305306221819d+00, &
   0.6877299149965d-08, 0.2834987233449d+01, 0.9561746721300d-02, &
   0.8311234227190d-08, 0.2202951835758d+01, 0.3801276407308d+00, &
   0.6599472832414d-08, 0.4478581462618d+01, 0.1063314406849d+01, &
   0.6160491096549d-08, 0.5145858696411d+01, 0.1368660381889d+01, &
   0.6164772043891d-08, 0.3762976697911d+00, 0.4234171675140d+00, &
   0.6363248684450d-08, 0.3162246718685d+01, 0.1253008786510d-01, &
   0.6448587520999d-08, 0.3442693302119d+01, 0.5287268506303d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 141, NS0X) /&
   0.6431662283977d-08, 0.8977549136606d+00, 0.5306550935933d+00, &
   0.6351223158474d-08, 0.4306447410369d+01, 0.5217580628120d+02, &
   0.5476721393451d-08, 0.3888529177855d+01, 0.2221856701002d+01, &
   0.5341772572619d-08, 0.2655560662512d+01, 0.7466759693650d-01, &
   0.5337055758302d-08, 0.5164990735946d+01, 0.7489573444450d-01, &
   0.5373120816787d-08, 0.6041214553456d+01, 0.1274714967946d+00, &
   0.5392351705426d-08, 0.9177763485932d+00, 0.1055449481598d+01, &
   0.6688495850205d-08, 0.3089608126937d+01, 0.2213766559277d+00, &
   0.5072003660362d-08, 0.4311316541553d+01, 0.2132517061319d+00, &
   0.9698639532817d-09, 0.1074134313634d+01, 0.7826370942180d+02, &
   0.5070726650455d-08, 0.5790675464444d+00, 0.2133464534247d+00 /
   !  SS-to-Sun, T^1, X
   DATA ((S1(I, J, 1), I = 1, 3), J = 1, 10) /&
   -0.1296310361520d-07, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.8975769009438d-08, 0.1128891609250d+01, 0.4265981595566d+00, &
   0.7771113441307d-08, 0.2706039877077d+01, 0.2061856251104d+00, &
   0.7538303866642d-08, 0.2191281289498d+01, 0.2204125344462d+00, &
   0.6061384579336d-08, 0.3248167319958d+01, 0.1059381944224d+01, &
   0.5726994235594d-08, 0.5569981398610d+01, 0.5225775174439d+00, &
   0.5616492836424d-08, 0.5057386614909d+01, 0.5368044267797d+00, &
   0.1010881584769d-08, 0.3473577116095d+01, 0.7113454667900d-02, &
   0.7259606157626d-09, 0.3651858593665d+00, 0.6398972393349d+00, &
   0.8755095026935d-09, 0.1662835408338d+01, 0.4194847048887d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 11, 20) /&
   0.5370491182812d-09, 0.1327673878077d+01, 0.4337116142245d+00, &
   0.5743773887665d-09, 0.4250200846687d+01, 0.2132990797783d+00, &
   0.4408103140300d-09, 0.3598752574277d+01, 0.1589072916335d+01, &
   0.3101892374445d-09, 0.4887822983319d+01, 0.1052268489556d+01, &
   0.3209453713578d-09, 0.9702272295114d+00, 0.5296909721118d+00, &
   0.3017228286064d-09, 0.5484462275949d+01, 0.1066495398892d+01, &
   0.3200700038601d-09, 0.2846613338643d+01, 0.1495633313810d+00, &
   0.2137637279911d-09, 0.5692163292729d+00, 0.3163918923335d+00, &
   0.1899686386727d-09, 0.2061077157189d+01, 0.2275259891141d+00, &
   0.1401994545308d-09, 0.4177771136967d+01, 0.1102062672231d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 21, 30) /&
   0.1578057810499d-09, 0.5782460597335d+01, 0.7626583626240d-01, &
   0.1237713253351d-09, 0.5705900866881d+01, 0.5154640627760d+00, &
   0.1313076837395d-09, 0.5163438179576d+01, 0.3664874755930d-01, &
   0.1184963304860d-09, 0.3054804427242d+01, 0.6327837846670d+00, &
   0.1238130878565d-09, 0.2317292575962d+01, 0.3961708870310d-01, &
   0.1015959527736d-09, 0.2194643645526d+01, 0.7329749511860d-01, &
   0.9017954423714d-10, 0.2868603545435d+01, 0.1990721704425d+00, &
   0.8668024955603d-10, 0.4923849675082d+01, 0.5439178814476d+00, &
   0.7756083930103d-10, 0.3014334135200d+01, 0.9491756770005d+00, &
   0.7536503401741d-10, 0.2704886279769d+01, 0.1030928125552d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 31, NS1X) /&
   0.5483308679332d-10, 0.6010983673799d+01, 0.8531963191132d+00, &
   0.5184339620428d-10, 0.1952704573291d+01, 0.2093666171530d+00, &
   0.5108658712030d-10, 0.2958575786649d+01, 0.2172315424036d+00, &
   0.5019424524650d-10, 0.1736317621318d+01, 0.2164800718209d+00, &
   0.4909312625978d-10, 0.3167216416257d+01, 0.2101180877357d+00, &
   0.4456638901107d-10, 0.7697579923471d+00, 0.3235053470014d+00, &
   0.4227030350925d-10, 0.3490910137928d+01, 0.6373574839730d-01, &
   0.4095456040093d-10, 0.5178888984491d+00, 0.6470106940028d+00, &
   0.4990537041422d-10, 0.3323887668974d+01, 0.1422690933580d-01, &
   0.4321170010845d-10, 0.4288484987118d+01, 0.7358765972222d+00 /

   !  SS-to-Sun, T^2, X
   DATA ((S2(I, J, 1), I = 1, 3), J = 1, NS2X) /&
   0.1603551636587d-11, 0.4404109410481d+01, 0.2061856251104d+00, &
   0.1556935889384d-11, 0.4818040873603d+00, 0.2204125344462d+00, &
   0.1182594414915d-11, 0.9935762734472d+00, 0.5225775174439d+00, &
   0.1158794583180d-11, 0.3353180966450d+01, 0.5368044267797d+00, &
   0.9597358943932d-12, 0.5567045358298d+01, 0.2132990797783d+00, &
   0.6511516579605d-12, 0.5630872420788d+01, 0.4265981595566d+00, &
   0.7419792747688d-12, 0.2156188581957d+01, 0.5296909721118d+00, &
   0.3951972655848d-12, 0.1981022541805d+01, 0.1059381944224d+01, &
   0.4478223877045d-12, 0.0000000000000d+00, 0.0000000000000d+00 /

   !  SS-to-Sun, T^0, Y
   DATA ((S0(I, J, 2), I = 1, 3), J = 1, 10) /&
   0.4955392320126d-02, 0.2170467313679d+01, 0.5296909721118d+00, &
   0.2722325167392d-02, 0.2444433682196d+01, 0.2132990797783d+00, &
   0.1546579925346d-02, 0.5992779281546d+00, 0.3813291813120d-01, &
   0.8363140252966d-03, 0.7687356310801d+00, 0.7478166569050d-01, &
   0.3385792683603d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.1201192221613d-03, 0.2520035601514d+01, 0.1059381944224d+01, &
   0.7587125720554d-04, 0.1669954006449d+01, 0.4265981595566d+00, &
   0.1964155361250d-04, 0.5707743963343d+01, 0.2061856251104d+00, &
   0.1891900364909d-04, 0.2320960679937d+01, 0.2204125344462d+00, &
   0.1937373433356d-04, 0.3226940689555d+01, 0.1495633313810d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 11, 20) /&
   0.1437139941351d-04, 0.2301626908096d+01, 0.5225775174439d+00, &
   0.1406267683099d-04, 0.5188579265542d+01, 0.5368044267797d+00, &
   0.1178703080346d-04, 0.5489483248476d+01, 0.7626583626240d-01, &
   0.8079835186041d-05, 0.1683751835264d+01, 0.3664874755930d-01, &
   0.7623253594652d-05, 0.2656400462961d+01, 0.3961708870310d-01, &
   0.6248667483971d-05, 0.4992775362055d+01, 0.7329749511860d-01, &
   0.4366353695038d-05, 0.2869706279678d+01, 0.1589072916335d+01, &
   0.3829101568895d-05, 0.3572131359950d+01, 0.7113454667900d-02, &
   0.3175733773908d-05, 0.4535372530045d+01, 0.4194847048887d+00, &
   0.3092437902159d-05, 0.9230153317909d+00, 0.6398972393349d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 21, 30) /&
   0.2874168812154d-05, 0.3363143761101d+01, 0.1102062672231d+00, &
   0.3040119321826d-05, 0.3324250895675d+01, 0.6283075850446d+01, &
   0.2699723308006d-05, 0.2917882441928d+00, 0.1030928125552d+00, &
   0.2134832683534d-05, 0.4220997202487d+01, 0.3163918923335d+00, &
   0.1770412139433d-05, 0.4747318496462d+01, 0.1021328554739d+02, &
   0.1377264209373d-05, 0.4305058462401d+00, 0.1484170571900d-02, &
   0.1127814538960d-05, 0.8538177240740d+00, 0.6327837846670d+00, &
   0.1055608090130d-05, 0.1551800742580d+01, 0.4337116142245d+00, &
   0.9802673861420d-06, 0.1459646735377d+01, 0.1052268489556d+01, &
   0.1090329461951d-05, 0.1587351228711d+01, 0.1162474756779d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 31, 40) /&
   0.6959590025090d-06, 0.5534442628766d+01, 0.1066495398892d+01, &
   0.5664914529542d-06, 0.6030673003297d+01, 0.9491756770005d+00, &
   0.6607787763599d-06, 0.4989507233927d+01, 0.8460828644453d+00, &
   0.6269725742838d-06, 0.4222951804572d+01, 0.1480791608091d+00, &
   0.6301889697863d-06, 0.5444316669126d+01, 0.2243449970715d+00, &
   0.4891042662861d-06, 0.1490552839784d+01, 0.3340612434717d+01, &
   0.3457083123290d-06, 0.3030475486049d+01, 0.3516457698740d-01, &
   0.3032559967314d-06, 0.2652038793632d+01, 0.1104591729320d-01, &
   0.2841133988903d-06, 0.1276744786829d+01, 0.4110125927500d-01, &
   0.2855564444432d-06, 0.2143368674733d+01, 0.1510475019529d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 41, 50) /&
   0.2765157135038d-06, 0.5444186109077d+01, 0.6373574839730d-01, &
   0.2382312465034d-06, 0.2190521137593d+01, 0.2275259891141d+00, &
   0.2808060365077d-06, 0.5735195064841d+01, 0.2535050500000d-01, &
   0.2332175234405d-06, 0.9481985524859d-01, 0.7181332454670d-01, &
   0.2322488199659d-06, 0.5180499361533d+01, 0.8582758298370d-01, &
   0.1881850258423d-06, 0.3219788273885d+01, 0.2118763888447d+01, &
   0.2196111392808d-06, 0.2366941159761d+01, 0.2968341143800d-02, &
   0.2183810335519d-06, 0.4825445110915d+01, 0.7775000683430d-01, &
   0.2002733093326d-06, 0.2457148995307d+01, 0.2093666171530d+00, &
   0.1967111767229d-06, 0.5586291545459d+01, 0.2172315424036d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 51, 60) /&
   0.1568473250543d-06, 0.3708003123320d+01, 0.7429900518901d+00, &
   0.1852528314300d-06, 0.4310638151560d+01, 0.2022531624851d+00, &
   0.1832111226447d-06, 0.1494665322656d+01, 0.3235053470014d+00, &
   0.1746805502310d-06, 0.1451378500784d+01, 0.1385174140878d+00, &
   0.1555730966650d-06, 0.1068040418198d+01, 0.7358765972222d+00, &
   0.1554883462559d-06, 0.2442579035461d+01, 0.5154640627760d+00, &
   0.1638380568746d-06, 0.2597913420625d+00, 0.8531963191132d+00, &
   0.1159938593640d-06, 0.5834512021280d+01, 0.1990721704425d+00, &
   0.1083427965695d-06, 0.5054033177950d+01, 0.5439178814476d+00, &
   0.1156480369431d-06, 0.5325677432457d+01, 0.5257585094865d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 61, 70) /&
   0.1141308860095d-06, 0.2153403923857d+01, 0.5336234347371d+00, &
   0.7913146470946d-07, 0.8642846847027d+00, 0.1478866649112d+01, &
   0.7439752463733d-07, 0.1970628496213d+01, 0.2164800718209d+00, &
   0.7280277104079d-07, 0.6073307250609d+01, 0.2101180877357d+00, &
   0.8319567719136d-07, 0.1954371928334d+01, 0.1692165728891d+01, &
   0.7137705549290d-07, 0.8904989440909d+00, 0.4155522422634d+00, &
   0.6900825396225d-07, 0.2825717714977d+01, 0.1173197218910d+00, &
   0.7245757216635d-07, 0.2481677513331d+01, 0.1265567569334d+01, &
   0.6961165696255d-07, 0.1292955312978d+01, 0.9562891316684d+00, &
   0.7571804456890d-07, 0.3427517575069d+01, 0.1422690933580d-01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 71, 80) /&
   0.6605425721904d-07, 0.8052192701492d+00, 0.6470106940028d+00, &
   0.7375477357248d-07, 0.1705076390088d+01, 0.1581959461667d+01, &
   0.7041664951470d-07, 0.4848356967891d+00, 0.9597935788730d-01, &
   0.6322199535763d-07, 0.3878069473909d+01, 0.7084920306520d-01, &
   0.5244380279191d-07, 0.2645560544125d+01, 0.5265099800692d+00, &
   0.5143125704988d-07, 0.4834486101370d+01, 0.5328719641544d+00, &
   0.5871866319373d-07, 0.7981472548900d+00, 0.7871412831580d-01, &
   0.6300822573871d-07, 0.5979398788281d+01, 0.2608790314060d+02, &
   0.6062154271548d-07, 0.4108655402756d+01, 0.1114304132498d+00, &
   0.4361912339976d-07, 0.5322624319280d+01, 0.1375773836557d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 81, 90) /&
   0.4417005920067d-07, 0.6240817359284d+01, 0.2770348281756d+00, &
   0.4686806749936d-07, 0.3214977301156d+01, 0.1143987543936d+00, &
   0.3758892132305d-07, 0.5879809634765d+01, 0.1596186371003d+01, &
   0.5151351332319d-07, 0.2893377688007d+00, 0.2228608264996d+00, &
   0.4554683578572d-07, 0.5475427144122d+01, 0.1465949902372d+00, &
   0.3442381385338d-07, 0.5992034796640d+01, 0.5070101000000d-01, &
   0.2831093954933d-07, 0.5367350273914d+01, 0.3092784376656d+00, &
   0.3756267090084d-07, 0.5758171285420d+01, 0.4903339079539d+00, &
   0.2816374679892d-07, 0.1863718700923d+01, 0.2991266627620d+00, &
   0.3419307025569d-07, 0.9524347534130d+00, 0.3518164938661d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 91, 100) /&
   0.2904250494239d-07, 0.5304471615602d+01, 0.1099462426779d+00, &
   0.2471734511206d-07, 0.1297069793530d+01, 0.6256703299991d+00, &
   0.2539620831872d-07, 0.3281126083375d+01, 0.1256615170089d+02, &
   0.2281017868007d-07, 0.1829122133165d+01, 0.6681224869435d+01, &
   0.2275319473335d-07, 0.5797198160181d+01, 0.3932462625300d-02, &
   0.2547755368442d-07, 0.4752697708330d+01, 0.1169588211447d+01, &
   0.2285979669317d-07, 0.1223205292886d+01, 0.1045155034888d+01, &
   0.1913386560994d-07, 0.1757532993389d+01, 0.1155361302111d+01, &
   0.1809020525147d-07, 0.4246116108791d+01, 0.3368040641550d-01, &
   0.1649213300201d-07, 0.1445162890627d+01, 0.4408250688924d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 101, 110) /&
   0.1834972793932d-07, 0.1126917567225d+01, 0.4452511715700d-02, &
   0.1439550648138d-07, 0.6160756834764d+01, 0.9420622223326d+00, &
   0.1487645457041d-07, 0.4358761931792d+01, 0.4123712502208d+00, &
   0.1731729516660d-07, 0.6134456753344d+01, 0.2108507877249d+00, &
   0.1717747163567d-07, 0.1898186084455d+01, 0.2157473718317d+00, &
   0.1418190430374d-07, 0.4180286741266d+01, 0.6521991896920d-01, &
   0.1404844134873d-07, 0.7654053565412d-01, 0.4258542984690d-01, &
   0.1409842846538d-07, 0.4418612420312d+01, 0.2258291676434d+00, &
   0.1090948346291d-07, 0.1260615686131d+01, 0.4226656969313d+00, &
   0.1357577323612d-07, 0.3558248818690d+01, 0.7923417740620d-01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 111, 120) /&
   0.1018154061960d-07, 0.5676087241256d+01, 0.1456308687557d+00, &
   0.1412073972109d-07, 0.8394392632422d+00, 0.1525316725248d+00, &
   0.1030938326496d-07, 0.1653593274064d+01, 0.1795258541446d+01, &
   0.1180081567104d-07, 0.1285802592036d+01, 0.7032915397480d-01, &
   0.9708510575650d-08, 0.7631889488106d+00, 0.8434341241180d-01, &
   0.9637689663447d-08, 0.4630642649176d+01, 0.1272681024002d+01, &
   0.1068910429389d-07, 0.5294934032165d+01, 0.2123349582968d+00, &
   0.1063716179336d-07, 0.2736266800832d+01, 0.2142632012598d+00, &
   0.1234858713814d-07, 0.1302891146570d+01, 0.1847279083684d+00, &
   0.8912631189738d-08, 0.3570415993621d+01, 0.2648454860559d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 121, 130) /&
   0.1036378285534d-07, 0.4236693440949d+01, 0.1370332435159d+00, &
   0.9667798501561d-08, 0.2960768892398d+01, 0.4376440768498d+00, &
   0.8108314201902d-08, 0.6987781646841d+00, 0.2880807454688d+00, &
   0.7648364324628d-08, 0.2499017863863d+01, 0.2037373330570d+00, &
   0.7286136828406d-08, 0.3787426951665d+01, 0.1129145838217d+00, &
   0.9448237743913d-08, 0.2694354332983d+01, 0.5272426800584d+00, &
   0.9374276106428d-08, 0.4787121277064d+01, 0.5321392641652d+00, &
   0.7100226287462d-08, 0.3530238792101d+00, 0.6288513220417d+00, &
   0.9253056659571d-08, 0.1399478925664d+01, 0.1606092486742d+00, &
   0.6636432145504d-08, 0.3479575438447d+01, 0.1368660381889d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 131, 140) /&
   0.6469975312932d-08, 0.1383669964800d+01, 0.2008557621224d+01, &
   0.7335849729765d-08, 0.1243698166898d+01, 0.9561746721300d-02, &
   0.8743421205855d-08, 0.3776164289301d+01, 0.3801276407308d+00, &
   0.5993635744494d-08, 0.5627122113596d+01, 0.2042657109477d+02, &
   0.5981008479693d-08, 0.1674336636752d+01, 0.2111650433779d+01, &
   0.6188535145838d-08, 0.5214925208672d+01, 0.4305306221819d+00, &
   0.6596074017566d-08, 0.2907653268124d+01, 0.1063314406849d+01, &
   0.6630815126226d-08, 0.2127643669658d+01, 0.8389694097774d+00, &
   0.6156772830040d-08, 0.5082160803295d+01, 0.4234171675140d+00, &
   0.6446960563014d-08, 0.1872100916905d+01, 0.5287268506303d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 141, NS0Y) /&
   0.6429324424668d-08, 0.5610276103577d+01, 0.5306550935933d+00, &
   0.6302232396465d-08, 0.1592152049607d+01, 0.1253008786510d-01, &
   0.6399244436159d-08, 0.2746214421532d+01, 0.5217580628120d+02, &
   0.5474965172558d-08, 0.2317666374383d+01, 0.2221856701002d+01, &
   0.5339293190692d-08, 0.1084724961156d+01, 0.7466759693650d-01, &
   0.5334733683389d-08, 0.3594106067745d+01, 0.7489573444450d-01, &
   0.5392665782110d-08, 0.5630254365606d+01, 0.1055449481598d+01, &
   0.6682075673789d-08, 0.1518480041732d+01, 0.2213766559277d+00, &
   0.5079130495960d-08, 0.2739765115711d+01, 0.2132517061319d+00, &
   0.9757679038404d-09, 0.5796846023126d+01, 0.7826370942180d+02, &
   0.5077759793261d-08, 0.5290711290094d+01, 0.2133464534247d+00 /
   !  SS-to-Sun, T^1, Y
   DATA ((S1(I, J, 2), I = 1, 3), J = 1, 10) /&
   0.8989047573576d-08, 0.5840593672122d+01, 0.4265981595566d+00, &
   0.7815938401048d-08, 0.1129664707133d+01, 0.2061856251104d+00, &
   0.7550926713280d-08, 0.6196589104845d+00, 0.2204125344462d+00, &
   0.6056556925895d-08, 0.1677494667846d+01, 0.1059381944224d+01, &
   0.5734142698204d-08, 0.4000920852962d+01, 0.5225775174439d+00, &
   0.5614341822459d-08, 0.3486722577328d+01, 0.5368044267797d+00, &
   0.1028678147656d-08, 0.1877141024787d+01, 0.7113454667900d-02, &
   0.7270792075266d-09, 0.5077167301739d+01, 0.6398972393349d+00, &
   0.8734141726040d-09, 0.9069550282609d-01, 0.4194847048887d+00, &
   0.5377371402113d-09, 0.6039381844671d+01, 0.4337116142245d+00 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 11, 20) /&
   0.4729719431571d-09, 0.2153086311760d+01, 0.2132990797783d+00, &
   0.4458052820973d-09, 0.5059830025565d+01, 0.5296909721118d+00, &
   0.4406855467908d-09, 0.2027971692630d+01, 0.1589072916335d+01, &
   0.3101659310977d-09, 0.3317677981860d+01, 0.1052268489556d+01, &
   0.3016749232545d-09, 0.3913703482532d+01, 0.1066495398892d+01, &
   0.3198541352656d-09, 0.1275513098525d+01, 0.1495633313810d+00, &
   0.2142065389871d-09, 0.5301351614597d+01, 0.3163918923335d+00, &
   0.1902615247592d-09, 0.4894943352736d+00, 0.2275259891141d+00, &
   0.1613410990871d-09, 0.2449891130437d+01, 0.1102062672231d+00, &
   0.1576992165097d-09, 0.4211421447633d+01, 0.7626583626240d-01 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 21, 30) /&
   0.1241637259894d-09, 0.4140803368133d+01, 0.5154640627760d+00, &
   0.1313974830355d-09, 0.3591920305503d+01, 0.3664874755930d-01, &
   0.1181697118258d-09, 0.1506314382788d+01, 0.6327837846670d+00, &
   0.1238239742779d-09, 0.7461405378404d+00, 0.3961708870310d-01, &
   0.1010107068241d-09, 0.6271010795475d+00, 0.7329749511860d-01, &
   0.9226316616509d-10, 0.1259158839583d+01, 0.1990721704425d+00, &
   0.8664946419555d-10, 0.3353244696934d+01, 0.5439178814476d+00, &
   0.7757230468978d-10, 0.1447677295196d+01, 0.9491756770005d+00, &
   0.7693168628139d-10, 0.1120509896721d+01, 0.1030928125552d+00, &
   0.5487897454612d-10, 0.4439380426795d+01, 0.8531963191132d+00 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 31, NS1Y) /&
   0.5196118677218d-10, 0.3788856619137d+00, 0.2093666171530d+00, &
   0.5110853339935d-10, 0.1386879372016d+01, 0.2172315424036d+00, &
   0.5027804534813d-10, 0.1647881805466d+00, 0.2164800718209d+00, &
   0.4922485922674d-10, 0.1594315079862d+01, 0.2101180877357d+00, &
   0.6155599524400d-10, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.4447147832161d-10, 0.5480720918976d+01, 0.3235053470014d+00, &
   0.4144691276422d-10, 0.1931371033660d+01, 0.6373574839730d-01, &
   0.4099950625452d-10, 0.5229611294335d+01, 0.6470106940028d+00, &
   0.5060541682953d-10, 0.1731112486298d+01, 0.1422690933580d-01, &
   0.4293615946300d-10, 0.2714571038925d+01, 0.7358765972222d+00 /

   !  SS-to-Sun, T^2, Y
   DATA ((S2(I, J, 2), I = 1, 3), J = 1, NS2Y) /&
   0.1609114495091d-11, 0.2831096993481d+01, 0.2061856251104d+00, &
   0.1560330784946d-11, 0.5193058213906d+01, 0.2204125344462d+00, &
   0.1183535479202d-11, 0.5707003443890d+01, 0.5225775174439d+00, &
   0.1158183066182d-11, 0.1782400404928d+01, 0.5368044267797d+00, &
   0.1032868027407d-11, 0.4036925452011d+01, 0.2132990797783d+00, &
   0.6540142847741d-12, 0.4058241056717d+01, 0.4265981595566d+00, &
   0.7305236491596d-12, 0.6175401942957d+00, 0.5296909721118d+00, &
   -0.5580725052968d-12, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3946122651015d-12, 0.4108265279171d+00, 0.1059381944224d+01 /

   !  SS-to-Sun, T^0, Z
   DATA ((S0(I, J, 3), I = 1, 3), J = 1, 10) /&
   0.1181255122986d-03, 0.4607918989164d+00, 0.2132990797783d+00, &
   0.1127777651095d-03, 0.4169146331296d+00, 0.5296909721118d+00, &
   0.4777754401806d-04, 0.4582657007130d+01, 0.3813291813120d-01, &
   0.1129354285772d-04, 0.5758735142480d+01, 0.7478166569050d-01, &
   -0.1149543637123d-04, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3298730512306d-05, 0.5978801994625d+01, 0.4265981595566d+00, &
   0.2733376706079d-05, 0.7665413691040d+00, 0.1059381944224d+01, &
   0.9426389657270d-06, 0.3710201265838d+01, 0.2061856251104d+00, &
   0.8187517749552d-06, 0.3390675605802d+00, 0.2204125344462d+00, &
   0.4080447871819d-06, 0.4552296640088d+00, 0.5225775174439d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 11, 20) /&
   0.3169973017028d-06, 0.3445455899321d+01, 0.5368044267797d+00, &
   0.2438098615549d-06, 0.5664675150648d+01, 0.3664874755930d-01, &
   0.2601897517235d-06, 0.1931894095697d+01, 0.1495633313810d+00, &
   0.2314558080079d-06, 0.3666319115574d+00, 0.3961708870310d-01, &
   0.1962549548002d-06, 0.3167411699020d+01, 0.7626583626240d-01, &
   0.2180518287925d-06, 0.1544420746580d+01, 0.7113454667900d-02, &
   0.1451382442868d-06, 0.1583756740070d+01, 0.1102062672231d+00, &
   0.1358439007389d-06, 0.5239941758280d+01, 0.6398972393349d+00, &
   0.1050585898028d-06, 0.2266958352859d+01, 0.3163918923335d+00, &
   0.1050029870186d-06, 0.2711495250354d+01, 0.4194847048887d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 21, 30) /&
   0.9934920679800d-07, 0.1116208151396d+01, 0.1589072916335d+01, &
   0.1048395331560d-06, 0.3408619600206d+01, 0.1021328554739d+02, &
   0.8370147196668d-07, 0.3810459401087d+01, 0.2535050500000d-01, &
   0.7989856510998d-07, 0.3769910473647d+01, 0.7329749511860d-01, &
   0.5441221655233d-07, 0.2416994903374d+01, 0.1030928125552d+00, &
   0.4610812906784d-07, 0.5858503336994d+01, 0.4337116142245d+00, &
   0.3923022803444d-07, 0.3354170010125d+00, 0.1484170571900d-02, &
   0.2610725582128d-07, 0.5410600646324d+01, 0.6327837846670d+00, &
   0.2455279767721d-07, 0.6120216681403d+01, 0.1162474756779d+01, &
   0.2375530706525d-07, 0.6055443426143d+01, 0.1052268489556d+01 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 31, 40) /&
   0.1782967577553d-07, 0.3146108708004d+01, 0.8460828644453d+00, &
   0.1581687095238d-07, 0.6255496089819d+00, 0.3340612434717d+01, &
   0.1594657672461d-07, 0.3782604300261d+01, 0.1066495398892d+01, &
   0.1563448615040d-07, 0.1997775733196d+01, 0.2022531624851d+00, &
   0.1463624258525d-07, 0.1736316792088d+00, 0.3516457698740d-01, &
   0.1331585056673d-07, 0.4331941830747d+01, 0.9491756770005d+00, &
   0.1130634557637d-07, 0.6152017751825d+01, 0.2968341143800d-02, &
   0.1028949607145d-07, 0.2101792614637d+00, 0.2275259891141d+00, &
   0.1024074971618d-07, 0.4071833211074d+01, 0.5070101000000d-01, &
   0.8826956060303d-08, 0.4861633688145d+00, 0.2093666171530d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 41, NS0Z) /&
   0.8572230171541d-08, 0.5268190724302d+01, 0.4110125927500d-01, &
   0.7649332643544d-08, 0.5134543417106d+01, 0.2608790314060d+02, &
   0.8581673291033d-08, 0.2920218146681d+01, 0.1480791608091d+00, &
   0.8430589300938d-08, 0.3604576619108d+01, 0.2172315424036d+00, &
   0.7776165501012d-08, 0.3772942249792d+01, 0.6373574839730d-01, &
   0.8311070234408d-08, 0.6200412329888d+01, 0.3235053470014d+00, &
   0.6927365212582d-08, 0.4543353113437d+01, 0.8531963191132d+00, &
   0.6791574208598d-08, 0.2882188406238d+01, 0.7181332454670d-01, &
   0.5593100811839d-08, 0.1776646892780d+01, 0.7429900518901d+00, &
   0.7788777276590d-09, 0.1900569908215d+01, 0.5217580628120d+02, &
   0.4553381853021d-08, 0.3949617611240d+01, 0.7775000683430d-01 /

   !  SS-to-Sun, T^1, Z
   DATA ((S1(I, J, 3), I = 1, 3), J = 1, 10) /&
   0.5444220475678d-08, 0.1803825509310d+01, 0.2132990797783d+00, &
   0.3883412695596d-08, 0.4668616389392d+01, 0.5296909721118d+00, &
   0.1334341434551d-08, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3730001266883d-09, 0.5401405918943d+01, 0.2061856251104d+00, &
   0.2894929197956d-09, 0.4932415609852d+01, 0.2204125344462d+00, &
   0.2857950357701d-09, 0.3154625362131d+01, 0.7478166569050d-01, &
   0.2499226432292d-09, 0.3657486128988d+01, 0.4265981595566d+00, &
   0.1937705443593d-09, 0.5740434679002d+01, 0.1059381944224d+01, &
   0.1374894396320d-09, 0.1712857366891d+01, 0.5368044267797d+00, &
   0.1217248678408d-09, 0.2312090870932d+01, 0.5225775174439d+00 /
   DATA ((S1(I, J, 3), I = 1, 3), J = 11, NS1Z) /&
   0.7961052740870d-10, 0.5283368554163d+01, 0.3813291813120d-01, &
   0.4979225949689d-10, 0.4298290471860d+01, 0.4194847048887d+00, &
   0.4388552286597d-10, 0.6145515047406d+01, 0.7113454667900d-02, &
   0.2586835212560d-10, 0.3019448001809d+01, 0.6398972393349d+00 /

   !  SS-to-Sun, T^2, Z
   DATA ((S2(I, J, 3), I = 1, 3), J = 1, NS2Z) /&
   0.3749920358054d-12, 0.3230285558668d+01, 0.2132990797783d+00, &
   0.2735037220939d-12, 0.6154322683046d+01, 0.5296909721118d+00 /

   t = (xjd-xjd2000)/yearj
   t2= t * t
   do j = 1, 3 ! loop on x, y, z
      aa = 0d0
      bb = 0d0

      do k = 1, ns0(j)
         arg = s0(2, k, j) + s0(3, k, j) * t
         amp = s0(1, k, j)
         freq = s0(3, k, j)
         ampv = freq * amp

         if (abs(amp) > poseps) then
            aa = aa + amp * cos(arg)
         endif

         if (abs(ampv) > viteps) then
            bb = bb - ampv * sin(arg)
         endif

      enddo

      do k = 1, ns1(j)
         arg = s1(2, k, j) + s1(3, k, j) * t
         amp = s1(1, k, j)
         freq = s1(3, k, j)

         if (abs(100d0 * amp) > poseps) then !effect over 100 years has been considered
            aa = aa + amp * t * cos(arg)
         endif

         if ((abs(amp) > viteps).or.(abs(100d0 * amp * freq) > viteps)) then !effect over 100 years has been considered
            bb = bb + (amp * cos(arg) - freq * amp * t * sin(arg))
         endif

      enddo

      do k = 1, ns2(j)
         arg = s2(2, k, j) + s2(3, k, j) * t
         cc = cos(arg)
         amp = s2(1, k, j)
         aa = aa + amp * t2 * cc
         bb = bb + (2d0 * t * amp * cc - s2(3, k, j) * amp * t2 * sin(arg))
      enddo

      xyz(j) = aa ! position in au
      xyzp(j) = bb/yearj ! velocity in au/day
   enddo
   return
end subroutine baryc_m
!********************************************************
subroutine baryc_hm(xjd, xyz, xyzp)
   !********************************************************
   !
   !     Vector Barycenter to  Sun at the julian date xjd
   !     xyz(i)  in au,  xyzp in au/day
   !     ecliptic frame at J2000.
   !     Highest accuracy
   !
   !     Source ! from VSOP Series provided by P. Bretagnon June 2000
   !
   !     F. Mignard OCA/Cassiopee
   !
   !    *INPUT
   !     xjd     : julian day
   !
   !    *OUTPUT
   !     xyz     : Sun from SS barycenter in au in ecliptic frame J2000
   !     xyzp    : Velocity of the Sun in au/day
   !
   !     Precision : 3.0e-8  au      ~  5  km
   !                 6.0e-10 au/day  ~  1.0 mm/s
   !
   !********************************************************
   !

   implicit none

   integer, parameter :: ns0m = 151, ns1m = 40, ns2m = 9
   integer, parameter :: ns0x = 151, ns1x = 40, ns2x = 9
   integer, parameter :: ns0y = 151, ns1y = 40, ns2y = 9
   integer, parameter :: ns0z = 51, ns1z = 14, ns2z = 2

   real(kind = dp), intent(in) :: xjd
   real(kind = dp), dimension(3), intent(out) :: xyz, xyzp

   integer :: ns0(3) = (/ns0x, ns0y, ns0z/)
   integer :: ns1(3) = (/ns1x, ns1y, ns1z/)
   integer :: ns2(3) = (/ns2x, ns2y, ns2z/)

   real(kind = dp) :: s0(3, ns0m, 3), s1(3, ns1m, 3), s2(3, ns2m, 3)

   real(kind = dp) :: t, t2, aa, bb, arg, amp, cc
   integer :: i, j, k


   !SSB-to-Sun, T^0, X
   DATA ((S0(I, J, 1), I = 1, 3), J = 1, 10) /&
   0.4956757536410d-02, 0.3741073751789d+01, 0.5296909721118d+00, &
   0.2718490072522d-02, 0.4016011511425d+01, 0.2132990797783d+00, &
   0.1546493974344d-02, 0.2170528330642d+01, 0.3813291813120d-01, &
   0.8366855276341d-03, 0.2339614075294d+01, 0.7478166569050d-01, &
   0.2936777942117d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.1201317439469d-03, 0.4090736353305d+01, 0.1059381944224d+01, &
   0.7578550887230d-04, 0.3241518088140d+01, 0.4265981595566d+00, &
   0.1941787367773d-04, 0.1012202064330d+01, 0.2061856251104d+00, &
   0.1889227765991d-04, 0.3892520416440d+01, 0.2204125344462d+00, &
   0.1937896968613d-04, 0.4797779441161d+01, 0.1495633313810d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 11, 20) /&
   0.1434506110873d-04, 0.3868960697933d+01, 0.5225775174439d+00, &
   0.1406659911580d-04, 0.4759766557397d+00, 0.5368044267797d+00, &
   0.1179022300202d-04, 0.7774961520598d+00, 0.7626583626240d-01, &
   0.8085864460959d-05, 0.3254654471465d+01, 0.3664874755930d-01, &
   0.7622752967615d-05, 0.4227633103489d+01, 0.3961708870310d-01, &
   0.6209171139066d-05, 0.2791828325711d+00, 0.7329749511860d-01, &
   0.4366435633970d-05, 0.4440454875925d+01, 0.1589072916335d+01, &
   0.3792124889348d-05, 0.5156393842356d+01, 0.7113454667900d-02, &
   0.3154548963402d-05, 0.6157005730093d+01, 0.4194847048887d+00, &
   0.3088359882942d-05, 0.2494567553163d+01, 0.6398972393349d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 21, 30) /&
   0.2788440902136d-05, 0.4934318747989d+01, 0.1102062672231d+00, &
   0.3039928456376d-05, 0.4895077702640d+01, 0.6283075850446d+01, &
   0.2272258457679d-05, 0.5278394064764d+01, 0.1030928125552d+00, &
   0.2162007057957d-05, 0.5802978019099d+01, 0.3163918923335d+00, &
   0.1767632855737d-05, 0.3415346595193d-01, 0.1021328554739d+02, &
   0.1349413459362d-05, 0.2001643230755d+01, 0.1484170571900d-02, &
   0.1170141900476d-05, 0.2424750491620d+01, 0.6327837846670d+00, &
   0.1054355266820d-05, 0.3123311487576d+01, 0.4337116142245d+00, &
   0.9800822461610d-06, 0.3026258088130d+01, 0.1052268489556d+01, &
   0.1091203749931d-05, 0.3157811670347d+01, 0.1162474756779d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 31, 40) /&
   0.6960236715913d-06, 0.8219570542313d+00, 0.1066495398892d+01, &
   0.5689257296909d-06, 0.1323052375236d+01, 0.9491756770005d+00, &
   0.6613172135802d-06, 0.2765348881598d+00, 0.8460828644453d+00, &
   0.6277702517571d-06, 0.5794064466382d+01, 0.1480791608091d+00, &
   0.6304884066699d-06, 0.7323555380787d+00, 0.2243449970715d+00, &
   0.4897850467382d-06, 0.3062464235399d+01, 0.3340612434717d+01, &
   0.3759148598786d-06, 0.4588290469664d+01, 0.3516457698740d-01, &
   0.3110520548195d-06, 0.1374299536572d+01, 0.6373574839730d-01, &
   0.3064708359780d-06, 0.4222267485047d+01, 0.1104591729320d-01, &
   0.2856347168241d-06, 0.3714202944973d+01, 0.1510475019529d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 41, 50) /&
   0.2840945514288d-06, 0.2847972875882d+01, 0.4110125927500d-01, &
   0.2378951599405d-06, 0.3762072563388d+01, 0.2275259891141d+00, &
   0.2714229481417d-06, 0.1036049980031d+01, 0.2535050500000d-01, &
   0.2323551717307d-06, 0.4682388599076d+00, 0.8582758298370d-01, &
   0.1881790512219d-06, 0.4790565425418d+01, 0.2118763888447d+01, &
   0.2261353968371d-06, 0.1669144912212d+01, 0.7181332454670d-01, &
   0.2214546389848d-06, 0.3937717281614d+01, 0.2968341143800d-02, &
   0.2184915594933d-06, 0.1129169845099d+00, 0.7775000683430d-01, &
   0.2000164937936d-06, 0.4030009638488d+01, 0.2093666171530d+00, &
   0.1966105136719d-06, 0.8745955786834d+00, 0.2172315424036d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 51, 60) /&
   0.1904742332624d-06, 0.5919743598964d+01, 0.2022531624851d+00, &
   0.1657399705031d-06, 0.2549141484884d+01, 0.7358765972222d+00, &
   0.1574070533987d-06, 0.5277533020230d+01, 0.7429900518901d+00, &
   0.1832261651039d-06, 0.3064688127777d+01, 0.3235053470014d+00, &
   0.1733615346569d-06, 0.3011432799094d+01, 0.1385174140878d+00, &
   0.1549124014496d-06, 0.4005569132359d+01, 0.5154640627760d+00, &
   0.1637044713838d-06, 0.1831375966632d+01, 0.8531963191132d+00, &
   0.1123420082383d-06, 0.1180270407578d+01, 0.1990721704425d+00, &
   0.1083754165740d-06, 0.3414101320863d+00, 0.5439178814476d+00, &
   0.1156638012655d-06, 0.6130479452594d+00, 0.5257585094865d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 61, 70) /&
   0.1142548785134d-06, 0.3724761948846d+01, 0.5336234347371d+00, &
   0.7921463895965d-07, 0.2435425589361d+01, 0.1478866649112d+01, &
   0.7428600285231d-07, 0.3542144398753d+01, 0.2164800718209d+00, &
   0.8323211246747d-07, 0.3525058072354d+01, 0.1692165728891d+01, &
   0.7257595116312d-07, 0.1364299431982d+01, 0.2101180877357d+00, &
   0.7111185833236d-07, 0.2460478875808d+01, 0.4155522422634d+00, &
   0.6868090383716d-07, 0.4397327670704d+01, 0.1173197218910d+00, &
   0.7226419974175d-07, 0.4042647308905d+01, 0.1265567569334d+01, &
   0.6955642383177d-07, 0.2865047906085d+01, 0.9562891316684d+00, &
   0.7492139296331d-07, 0.5014278994215d+01, 0.1422690933580d-01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 71, 80) /&
   0.6598363128857d-07, 0.2376730020492d+01, 0.6470106940028d+00, &
   0.7381147293385d-07, 0.3272990384244d+01, 0.1581959461667d+01, &
   0.6402909624032d-07, 0.5302290955138d+01, 0.9597935788730d-01, &
   0.6237454263857d-07, 0.5444144425332d+01, 0.7084920306520d-01, &
   0.5241198544016d-07, 0.4215359579205d+01, 0.5265099800692d+00, &
   0.5144463853918d-07, 0.1218916689916d+00, 0.5328719641544d+00, &
   0.5868164772299d-07, 0.2369402002213d+01, 0.7871412831580d-01, &
   0.6233195669151d-07, 0.1254922242403d+01, 0.2608790314060d+02, &
   0.6068463791422d-07, 0.5679713760431d+01, 0.1114304132498d+00, &
   0.4359361135065d-07, 0.6097219641646d+00, 0.1375773836557d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 81, 90) /&
   0.4686510366826d-07, 0.4786231041431d+01, 0.1143987543936d+00, &
   0.3758977287225d-07, 0.1167368068139d+01, 0.1596186371003d+01, &
   0.4282051974778d-07, 0.1519471064319d+01, 0.2770348281756d+00, &
   0.5153765386113d-07, 0.1860532322984d+01, 0.2228608264996d+00, &
   0.4575129387188d-07, 0.7632857887158d+00, 0.1465949902372d+00, &
   0.3326844933286d-07, 0.1298219485285d+01, 0.5070101000000d-01, &
   0.3748617450984d-07, 0.1046510321062d+01, 0.4903339079539d+00, &
   0.2816756661499d-07, 0.3434522346190d+01, 0.2991266627620d+00, &
   0.3412750405039d-07, 0.2523766270318d+01, 0.3518164938661d+00, &
   0.2655796761776d-07, 0.2904422260194d+01, 0.6256703299991d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 91, 100) /&
   0.2963597929458d-07, 0.5923900431149d+00, 0.1099462426779d+00, &
   0.2539523734781d-07, 0.4851947722567d+01, 0.1256615170089d+02, &
   0.2283087914139d-07, 0.3400498595496d+01, 0.6681224869435d+01, &
   0.2321309799331d-07, 0.5789099148673d+01, 0.3368040641550d-01, &
   0.2549657649750d-07, 0.3991856479792d-01, 0.1169588211447d+01, &
   0.2290462303977d-07, 0.2788567577052d+01, 0.1045155034888d+01, &
   0.1945398522914d-07, 0.3290896998176d+01, 0.1155361302111d+01, &
   0.1849171512638d-07, 0.2698060129367d+01, 0.4452511715700d-02, &
   0.1647199834254d-07, 0.3016735644085d+01, 0.4408250688924d+00, &
   0.1529530765273d-07, 0.5573043116178d+01, 0.6521991896920d-01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 101, 110) /&
   0.1433199339978d-07, 0.1481192356147d+01, 0.9420622223326d+00, &
   0.1729134193602d-07, 0.1422817538933d+01, 0.2108507877249d+00, &
   0.1716463931346d-07, 0.3469468901855d+01, 0.2157473718317d+00, &
   0.1391206061378d-07, 0.6122436220547d+01, 0.4123712502208d+00, &
   0.1404746661924d-07, 0.1647765641936d+01, 0.4258542984690d-01, &
   0.1410452399455d-07, 0.5989729161964d+01, 0.2258291676434d+00, &
   0.1089828772168d-07, 0.2833705509371d+01, 0.4226656969313d+00, &
   0.1047374564948d-07, 0.5090690007331d+00, 0.3092784376656d+00, &
   0.1358279126532d-07, 0.5128990262836d+01, 0.7923417740620d-01, &
   0.1020456476148d-07, 0.9632772880808d+00, 0.1456308687557d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 111, 120) /&
   0.1033428735328d-07, 0.3223779318418d+01, 0.1795258541446d+01, &
   0.1412435841540d-07, 0.2410271572721d+01, 0.1525316725248d+00, &
   0.9722759371574d-08, 0.2333531395690d+01, 0.8434341241180d-01, &
   0.9657334084704d-08, 0.6199270974168d+01, 0.1272681024002d+01, &
   0.1083641148690d-07, 0.2864222292929d+01, 0.7032915397480d-01, &
   0.1067318403838d-07, 0.5833458866568d+00, 0.2123349582968d+00, &
   0.1062366201976d-07, 0.4307753989494d+01, 0.2142632012598d+00, &
   0.1236364149266d-07, 0.2873917870593d+01, 0.1847279083684d+00, &
   0.1092759489593d-07, 0.2959887266733d+01, 0.1370332435159d+00, &
   0.8912069362899d-08, 0.5141213702562d+01, 0.2648454860559d+01 /

   DATA ((S0(I, J, 1), I = 1, 3), J = 121, 130) /&
   0.9656467707970d-08, 0.4532182462323d+01, 0.4376440768498d+00, &
   0.8098386150135d-08, 0.2268906338379d+01, 0.2880807454688d+00, &
   0.7857714675000d-08, 0.4055544260745d+01, 0.2037373330570d+00, &
   0.7288455940646d-08, 0.5357901655142d+01, 0.1129145838217d+00, &
   0.9450595950552d-08, 0.4264926963939d+01, 0.5272426800584d+00, &
   0.9381718247537d-08, 0.7489366976576d-01, 0.5321392641652d+00, &
   0.7079052646038d-08, 0.1923311052874d+01, 0.6288513220417d+00, &
   0.9259004415344d-08, 0.2970256853438d+01, 0.1606092486742d+00, &
   0.8259801499742d-08, 0.3327056314697d+01, 0.8389694097774d+00, &
   0.6476334355779d-08, 0.2954925505727d+01, 0.2008557621224d+01 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 131, 140) /&
   0.5984021492007d-08, 0.9138753105829d+00, 0.2042657109477d+02, &
   0.5989546863181d-08, 0.3244464082031d+01, 0.2111650433779d+01, &
   0.6233108606023d-08, 0.4995232638403d+00, 0.4305306221819d+00, &
   0.6877299149965d-08, 0.2834987233449d+01, 0.9561746721300d-02, &
   0.8311234227190d-08, 0.2202951835758d+01, 0.3801276407308d+00, &
   0.6599472832414d-08, 0.4478581462618d+01, 0.1063314406849d+01, &
   0.6160491096549d-08, 0.5145858696411d+01, 0.1368660381889d+01, &
   0.6164772043891d-08, 0.3762976697911d+00, 0.4234171675140d+00, &
   0.6363248684450d-08, 0.3162246718685d+01, 0.1253008786510d-01, &
   0.6448587520999d-08, 0.3442693302119d+01, 0.5287268506303d+00 /
   DATA ((S0(I, J, 1), I = 1, 3), J = 141, NS0X) /&
   0.6431662283977d-08, 0.8977549136606d+00, 0.5306550935933d+00, &
   0.6351223158474d-08, 0.4306447410369d+01, 0.5217580628120d+02, &
   0.5476721393451d-08, 0.3888529177855d+01, 0.2221856701002d+01, &
   0.5341772572619d-08, 0.2655560662512d+01, 0.7466759693650d-01, &
   0.5337055758302d-08, 0.5164990735946d+01, 0.7489573444450d-01, &
   0.5373120816787d-08, 0.6041214553456d+01, 0.1274714967946d+00, &
   0.5392351705426d-08, 0.9177763485932d+00, 0.1055449481598d+01, &
   0.6688495850205d-08, 0.3089608126937d+01, 0.2213766559277d+00, &
   0.5072003660362d-08, 0.4311316541553d+01, 0.2132517061319d+00, &
   0.9698639532817d-09, 0.1074134313634d+01, 0.7826370942180d+02, &
   0.5070726650455d-08, 0.5790675464444d+00, 0.2133464534247d+00 /
   !  SS-to-Sun, T^1, X
   DATA ((S1(I, J, 1), I = 1, 3), J = 1, 10) /&
   -0.1296310361520d-07, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.8975769009438d-08, 0.1128891609250d+01, 0.4265981595566d+00, &
   0.7771113441307d-08, 0.2706039877077d+01, 0.2061856251104d+00, &
   0.7538303866642d-08, 0.2191281289498d+01, 0.2204125344462d+00, &
   0.6061384579336d-08, 0.3248167319958d+01, 0.1059381944224d+01, &
   0.5726994235594d-08, 0.5569981398610d+01, 0.5225775174439d+00, &
   0.5616492836424d-08, 0.5057386614909d+01, 0.5368044267797d+00, &
   0.1010881584769d-08, 0.3473577116095d+01, 0.7113454667900d-02, &
   0.7259606157626d-09, 0.3651858593665d+00, 0.6398972393349d+00, &
   0.8755095026935d-09, 0.1662835408338d+01, 0.4194847048887d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 11, 20) /&
   0.5370491182812d-09, 0.1327673878077d+01, 0.4337116142245d+00, &
   0.5743773887665d-09, 0.4250200846687d+01, 0.2132990797783d+00, &
   0.4408103140300d-09, 0.3598752574277d+01, 0.1589072916335d+01, &
   0.3101892374445d-09, 0.4887822983319d+01, 0.1052268489556d+01, &
   0.3209453713578d-09, 0.9702272295114d+00, 0.5296909721118d+00, &
   0.3017228286064d-09, 0.5484462275949d+01, 0.1066495398892d+01, &
   0.3200700038601d-09, 0.2846613338643d+01, 0.1495633313810d+00, &
   0.2137637279911d-09, 0.5692163292729d+00, 0.3163918923335d+00, &
   0.1899686386727d-09, 0.2061077157189d+01, 0.2275259891141d+00, &
   0.1401994545308d-09, 0.4177771136967d+01, 0.1102062672231d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 21, 30) /&
   0.1578057810499d-09, 0.5782460597335d+01, 0.7626583626240d-01, &
   0.1237713253351d-09, 0.5705900866881d+01, 0.5154640627760d+00, &
   0.1313076837395d-09, 0.5163438179576d+01, 0.3664874755930d-01, &
   0.1184963304860d-09, 0.3054804427242d+01, 0.6327837846670d+00, &
   0.1238130878565d-09, 0.2317292575962d+01, 0.3961708870310d-01, &
   0.1015959527736d-09, 0.2194643645526d+01, 0.7329749511860d-01, &
   0.9017954423714d-10, 0.2868603545435d+01, 0.1990721704425d+00, &
   0.8668024955603d-10, 0.4923849675082d+01, 0.5439178814476d+00, &
   0.7756083930103d-10, 0.3014334135200d+01, 0.9491756770005d+00, &
   0.7536503401741d-10, 0.2704886279769d+01, 0.1030928125552d+00 /
   DATA ((S1(I, J, 1), I = 1, 3), J = 31, NS1X) /&
   0.5483308679332d-10, 0.6010983673799d+01, 0.8531963191132d+00, &
   0.5184339620428d-10, 0.1952704573291d+01, 0.2093666171530d+00, &
   0.5108658712030d-10, 0.2958575786649d+01, 0.2172315424036d+00, &
   0.5019424524650d-10, 0.1736317621318d+01, 0.2164800718209d+00, &
   0.4909312625978d-10, 0.3167216416257d+01, 0.2101180877357d+00, &
   0.4456638901107d-10, 0.7697579923471d+00, 0.3235053470014d+00, &
   0.4227030350925d-10, 0.3490910137928d+01, 0.6373574839730d-01, &
   0.4095456040093d-10, 0.5178888984491d+00, 0.6470106940028d+00, &
   0.4990537041422d-10, 0.3323887668974d+01, 0.1422690933580d-01, &
   0.4321170010845d-10, 0.4288484987118d+01, 0.7358765972222d+00 /

   !  SS-to-Sun, T^2, X
   DATA ((S2(I, J, 1), I = 1, 3), J = 1, NS2X) /&
   0.1603551636587d-11, 0.4404109410481d+01, 0.2061856251104d+00, &
   0.1556935889384d-11, 0.4818040873603d+00, 0.2204125344462d+00, &
   0.1182594414915d-11, 0.9935762734472d+00, 0.5225775174439d+00, &
   0.1158794583180d-11, 0.3353180966450d+01, 0.5368044267797d+00, &
   0.9597358943932d-12, 0.5567045358298d+01, 0.2132990797783d+00, &
   0.6511516579605d-12, 0.5630872420788d+01, 0.4265981595566d+00, &
   0.7419792747688d-12, 0.2156188581957d+01, 0.5296909721118d+00, &
   0.3951972655848d-12, 0.1981022541805d+01, 0.1059381944224d+01, &
   0.4478223877045d-12, 0.0000000000000d+00, 0.0000000000000d+00 /

   !  SS-to-Sun, T^0, Y
   DATA ((S0(I, J, 2), I = 1, 3), J = 1, 10) /&
   0.4955392320126d-02, 0.2170467313679d+01, 0.5296909721118d+00, &
   0.2722325167392d-02, 0.2444433682196d+01, 0.2132990797783d+00, &
   0.1546579925346d-02, 0.5992779281546d+00, 0.3813291813120d-01, &
   0.8363140252966d-03, 0.7687356310801d+00, 0.7478166569050d-01, &
   0.3385792683603d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.1201192221613d-03, 0.2520035601514d+01, 0.1059381944224d+01, &
   0.7587125720554d-04, 0.1669954006449d+01, 0.4265981595566d+00, &
   0.1964155361250d-04, 0.5707743963343d+01, 0.2061856251104d+00, &
   0.1891900364909d-04, 0.2320960679937d+01, 0.2204125344462d+00, &
   0.1937373433356d-04, 0.3226940689555d+01, 0.1495633313810d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 11, 20) /&
   0.1437139941351d-04, 0.2301626908096d+01, 0.5225775174439d+00, &
   0.1406267683099d-04, 0.5188579265542d+01, 0.5368044267797d+00, &
   0.1178703080346d-04, 0.5489483248476d+01, 0.7626583626240d-01, &
   0.8079835186041d-05, 0.1683751835264d+01, 0.3664874755930d-01, &
   0.7623253594652d-05, 0.2656400462961d+01, 0.3961708870310d-01, &
   0.6248667483971d-05, 0.4992775362055d+01, 0.7329749511860d-01, &
   0.4366353695038d-05, 0.2869706279678d+01, 0.1589072916335d+01, &
   0.3829101568895d-05, 0.3572131359950d+01, 0.7113454667900d-02, &
   0.3175733773908d-05, 0.4535372530045d+01, 0.4194847048887d+00, &
   0.3092437902159d-05, 0.9230153317909d+00, 0.6398972393349d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 21, 30) /&
   0.2874168812154d-05, 0.3363143761101d+01, 0.1102062672231d+00, &
   0.3040119321826d-05, 0.3324250895675d+01, 0.6283075850446d+01, &
   0.2699723308006d-05, 0.2917882441928d+00, 0.1030928125552d+00, &
   0.2134832683534d-05, 0.4220997202487d+01, 0.3163918923335d+00, &
   0.1770412139433d-05, 0.4747318496462d+01, 0.1021328554739d+02, &
   0.1377264209373d-05, 0.4305058462401d+00, 0.1484170571900d-02, &
   0.1127814538960d-05, 0.8538177240740d+00, 0.6327837846670d+00, &
   0.1055608090130d-05, 0.1551800742580d+01, 0.4337116142245d+00, &
   0.9802673861420d-06, 0.1459646735377d+01, 0.1052268489556d+01, &
   0.1090329461951d-05, 0.1587351228711d+01, 0.1162474756779d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 31, 40) /&
   0.6959590025090d-06, 0.5534442628766d+01, 0.1066495398892d+01, &
   0.5664914529542d-06, 0.6030673003297d+01, 0.9491756770005d+00, &
   0.6607787763599d-06, 0.4989507233927d+01, 0.8460828644453d+00, &
   0.6269725742838d-06, 0.4222951804572d+01, 0.1480791608091d+00, &
   0.6301889697863d-06, 0.5444316669126d+01, 0.2243449970715d+00, &
   0.4891042662861d-06, 0.1490552839784d+01, 0.3340612434717d+01, &
   0.3457083123290d-06, 0.3030475486049d+01, 0.3516457698740d-01, &
   0.3032559967314d-06, 0.2652038793632d+01, 0.1104591729320d-01, &
   0.2841133988903d-06, 0.1276744786829d+01, 0.4110125927500d-01, &
   0.2855564444432d-06, 0.2143368674733d+01, 0.1510475019529d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 41, 50) /&
   0.2765157135038d-06, 0.5444186109077d+01, 0.6373574839730d-01, &
   0.2382312465034d-06, 0.2190521137593d+01, 0.2275259891141d+00, &
   0.2808060365077d-06, 0.5735195064841d+01, 0.2535050500000d-01, &
   0.2332175234405d-06, 0.9481985524859d-01, 0.7181332454670d-01, &
   0.2322488199659d-06, 0.5180499361533d+01, 0.8582758298370d-01, &
   0.1881850258423d-06, 0.3219788273885d+01, 0.2118763888447d+01, &
   0.2196111392808d-06, 0.2366941159761d+01, 0.2968341143800d-02, &
   0.2183810335519d-06, 0.4825445110915d+01, 0.7775000683430d-01, &
   0.2002733093326d-06, 0.2457148995307d+01, 0.2093666171530d+00, &
   0.1967111767229d-06, 0.5586291545459d+01, 0.2172315424036d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 51, 60) /&
   0.1568473250543d-06, 0.3708003123320d+01, 0.7429900518901d+00, &
   0.1852528314300d-06, 0.4310638151560d+01, 0.2022531624851d+00, &
   0.1832111226447d-06, 0.1494665322656d+01, 0.3235053470014d+00, &
   0.1746805502310d-06, 0.1451378500784d+01, 0.1385174140878d+00, &
   0.1555730966650d-06, 0.1068040418198d+01, 0.7358765972222d+00, &
   0.1554883462559d-06, 0.2442579035461d+01, 0.5154640627760d+00, &
   0.1638380568746d-06, 0.2597913420625d+00, 0.8531963191132d+00, &
   0.1159938593640d-06, 0.5834512021280d+01, 0.1990721704425d+00, &
   0.1083427965695d-06, 0.5054033177950d+01, 0.5439178814476d+00, &
   0.1156480369431d-06, 0.5325677432457d+01, 0.5257585094865d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 61, 70) /&
   0.1141308860095d-06, 0.2153403923857d+01, 0.5336234347371d+00, &
   0.7913146470946d-07, 0.8642846847027d+00, 0.1478866649112d+01, &
   0.7439752463733d-07, 0.1970628496213d+01, 0.2164800718209d+00, &
   0.7280277104079d-07, 0.6073307250609d+01, 0.2101180877357d+00, &
   0.8319567719136d-07, 0.1954371928334d+01, 0.1692165728891d+01, &
   0.7137705549290d-07, 0.8904989440909d+00, 0.4155522422634d+00, &
   0.6900825396225d-07, 0.2825717714977d+01, 0.1173197218910d+00, &
   0.7245757216635d-07, 0.2481677513331d+01, 0.1265567569334d+01, &
   0.6961165696255d-07, 0.1292955312978d+01, 0.9562891316684d+00, &
   0.7571804456890d-07, 0.3427517575069d+01, 0.1422690933580d-01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 71, 80) /&
   0.6605425721904d-07, 0.8052192701492d+00, 0.6470106940028d+00, &
   0.7375477357248d-07, 0.1705076390088d+01, 0.1581959461667d+01, &
   0.7041664951470d-07, 0.4848356967891d+00, 0.9597935788730d-01, &
   0.6322199535763d-07, 0.3878069473909d+01, 0.7084920306520d-01, &
   0.5244380279191d-07, 0.2645560544125d+01, 0.5265099800692d+00, &
   0.5143125704988d-07, 0.4834486101370d+01, 0.5328719641544d+00, &
   0.5871866319373d-07, 0.7981472548900d+00, 0.7871412831580d-01, &
   0.6300822573871d-07, 0.5979398788281d+01, 0.2608790314060d+02, &
   0.6062154271548d-07, 0.4108655402756d+01, 0.1114304132498d+00, &
   0.4361912339976d-07, 0.5322624319280d+01, 0.1375773836557d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 81, 90) /&
   0.4417005920067d-07, 0.6240817359284d+01, 0.2770348281756d+00, &
   0.4686806749936d-07, 0.3214977301156d+01, 0.1143987543936d+00, &
   0.3758892132305d-07, 0.5879809634765d+01, 0.1596186371003d+01, &
   0.5151351332319d-07, 0.2893377688007d+00, 0.2228608264996d+00, &
   0.4554683578572d-07, 0.5475427144122d+01, 0.1465949902372d+00, &
   0.3442381385338d-07, 0.5992034796640d+01, 0.5070101000000d-01, &
   0.2831093954933d-07, 0.5367350273914d+01, 0.3092784376656d+00, &
   0.3756267090084d-07, 0.5758171285420d+01, 0.4903339079539d+00, &
   0.2816374679892d-07, 0.1863718700923d+01, 0.2991266627620d+00, &
   0.3419307025569d-07, 0.9524347534130d+00, 0.3518164938661d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 91, 100) /&
   0.2904250494239d-07, 0.5304471615602d+01, 0.1099462426779d+00, &
   0.2471734511206d-07, 0.1297069793530d+01, 0.6256703299991d+00, &
   0.2539620831872d-07, 0.3281126083375d+01, 0.1256615170089d+02, &
   0.2281017868007d-07, 0.1829122133165d+01, 0.6681224869435d+01, &
   0.2275319473335d-07, 0.5797198160181d+01, 0.3932462625300d-02, &
   0.2547755368442d-07, 0.4752697708330d+01, 0.1169588211447d+01, &
   0.2285979669317d-07, 0.1223205292886d+01, 0.1045155034888d+01, &
   0.1913386560994d-07, 0.1757532993389d+01, 0.1155361302111d+01, &
   0.1809020525147d-07, 0.4246116108791d+01, 0.3368040641550d-01, &
   0.1649213300201d-07, 0.1445162890627d+01, 0.4408250688924d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 101, 110) /&
   0.1834972793932d-07, 0.1126917567225d+01, 0.4452511715700d-02, &
   0.1439550648138d-07, 0.6160756834764d+01, 0.9420622223326d+00, &
   0.1487645457041d-07, 0.4358761931792d+01, 0.4123712502208d+00, &
   0.1731729516660d-07, 0.6134456753344d+01, 0.2108507877249d+00, &
   0.1717747163567d-07, 0.1898186084455d+01, 0.2157473718317d+00, &
   0.1418190430374d-07, 0.4180286741266d+01, 0.6521991896920d-01, &
   0.1404844134873d-07, 0.7654053565412d-01, 0.4258542984690d-01, &
   0.1409842846538d-07, 0.4418612420312d+01, 0.2258291676434d+00, &
   0.1090948346291d-07, 0.1260615686131d+01, 0.4226656969313d+00, &
   0.1357577323612d-07, 0.3558248818690d+01, 0.7923417740620d-01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 111, 120) /&
   0.1018154061960d-07, 0.5676087241256d+01, 0.1456308687557d+00, &
   0.1412073972109d-07, 0.8394392632422d+00, 0.1525316725248d+00, &
   0.1030938326496d-07, 0.1653593274064d+01, 0.1795258541446d+01, &
   0.1180081567104d-07, 0.1285802592036d+01, 0.7032915397480d-01, &
   0.9708510575650d-08, 0.7631889488106d+00, 0.8434341241180d-01, &
   0.9637689663447d-08, 0.4630642649176d+01, 0.1272681024002d+01, &
   0.1068910429389d-07, 0.5294934032165d+01, 0.2123349582968d+00, &
   0.1063716179336d-07, 0.2736266800832d+01, 0.2142632012598d+00, &
   0.1234858713814d-07, 0.1302891146570d+01, 0.1847279083684d+00, &
   0.8912631189738d-08, 0.3570415993621d+01, 0.2648454860559d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 121, 130) /&
   0.1036378285534d-07, 0.4236693440949d+01, 0.1370332435159d+00, &
   0.9667798501561d-08, 0.2960768892398d+01, 0.4376440768498d+00, &
   0.8108314201902d-08, 0.6987781646841d+00, 0.2880807454688d+00, &
   0.7648364324628d-08, 0.2499017863863d+01, 0.2037373330570d+00, &
   0.7286136828406d-08, 0.3787426951665d+01, 0.1129145838217d+00, &
   0.9448237743913d-08, 0.2694354332983d+01, 0.5272426800584d+00, &
   0.9374276106428d-08, 0.4787121277064d+01, 0.5321392641652d+00, &
   0.7100226287462d-08, 0.3530238792101d+00, 0.6288513220417d+00, &
   0.9253056659571d-08, 0.1399478925664d+01, 0.1606092486742d+00, &
   0.6636432145504d-08, 0.3479575438447d+01, 0.1368660381889d+01 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 131, 140) /&
   0.6469975312932d-08, 0.1383669964800d+01, 0.2008557621224d+01, &
   0.7335849729765d-08, 0.1243698166898d+01, 0.9561746721300d-02, &
   0.8743421205855d-08, 0.3776164289301d+01, 0.3801276407308d+00, &
   0.5993635744494d-08, 0.5627122113596d+01, 0.2042657109477d+02, &
   0.5981008479693d-08, 0.1674336636752d+01, 0.2111650433779d+01, &
   0.6188535145838d-08, 0.5214925208672d+01, 0.4305306221819d+00, &
   0.6596074017566d-08, 0.2907653268124d+01, 0.1063314406849d+01, &
   0.6630815126226d-08, 0.2127643669658d+01, 0.8389694097774d+00, &
   0.6156772830040d-08, 0.5082160803295d+01, 0.4234171675140d+00, &
   0.6446960563014d-08, 0.1872100916905d+01, 0.5287268506303d+00 /
   DATA ((S0(I, J, 2), I = 1, 3), J = 141, NS0Y) /&
   0.6429324424668d-08, 0.5610276103577d+01, 0.5306550935933d+00, &
   0.6302232396465d-08, 0.1592152049607d+01, 0.1253008786510d-01, &
   0.6399244436159d-08, 0.2746214421532d+01, 0.5217580628120d+02, &
   0.5474965172558d-08, 0.2317666374383d+01, 0.2221856701002d+01, &
   0.5339293190692d-08, 0.1084724961156d+01, 0.7466759693650d-01, &
   0.5334733683389d-08, 0.3594106067745d+01, 0.7489573444450d-01, &
   0.5392665782110d-08, 0.5630254365606d+01, 0.1055449481598d+01, &
   0.6682075673789d-08, 0.1518480041732d+01, 0.2213766559277d+00, &
   0.5079130495960d-08, 0.2739765115711d+01, 0.2132517061319d+00, &
   0.9757679038404d-09, 0.5796846023126d+01, 0.7826370942180d+02, &
   0.5077759793261d-08, 0.5290711290094d+01, 0.2133464534247d+00 /
   !  SS-to-Sun, T^1, Y
   DATA ((S1(I, J, 2), I = 1, 3), J = 1, 10) /&
   0.8989047573576d-08, 0.5840593672122d+01, 0.4265981595566d+00, &
   0.7815938401048d-08, 0.1129664707133d+01, 0.2061856251104d+00, &
   0.7550926713280d-08, 0.6196589104845d+00, 0.2204125344462d+00, &
   0.6056556925895d-08, 0.1677494667846d+01, 0.1059381944224d+01, &
   0.5734142698204d-08, 0.4000920852962d+01, 0.5225775174439d+00, &
   0.5614341822459d-08, 0.3486722577328d+01, 0.5368044267797d+00, &
   0.1028678147656d-08, 0.1877141024787d+01, 0.7113454667900d-02, &
   0.7270792075266d-09, 0.5077167301739d+01, 0.6398972393349d+00, &
   0.8734141726040d-09, 0.9069550282609d-01, 0.4194847048887d+00, &
   0.5377371402113d-09, 0.6039381844671d+01, 0.4337116142245d+00 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 11, 20) /&
   0.4729719431571d-09, 0.2153086311760d+01, 0.2132990797783d+00, &
   0.4458052820973d-09, 0.5059830025565d+01, 0.5296909721118d+00, &
   0.4406855467908d-09, 0.2027971692630d+01, 0.1589072916335d+01, &
   0.3101659310977d-09, 0.3317677981860d+01, 0.1052268489556d+01, &
   0.3016749232545d-09, 0.3913703482532d+01, 0.1066495398892d+01, &
   0.3198541352656d-09, 0.1275513098525d+01, 0.1495633313810d+00, &
   0.2142065389871d-09, 0.5301351614597d+01, 0.3163918923335d+00, &
   0.1902615247592d-09, 0.4894943352736d+00, 0.2275259891141d+00, &
   0.1613410990871d-09, 0.2449891130437d+01, 0.1102062672231d+00, &
   0.1576992165097d-09, 0.4211421447633d+01, 0.7626583626240d-01 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 21, 30) /&
   0.1241637259894d-09, 0.4140803368133d+01, 0.5154640627760d+00, &
   0.1313974830355d-09, 0.3591920305503d+01, 0.3664874755930d-01, &
   0.1181697118258d-09, 0.1506314382788d+01, 0.6327837846670d+00, &
   0.1238239742779d-09, 0.7461405378404d+00, 0.3961708870310d-01, &
   0.1010107068241d-09, 0.6271010795475d+00, 0.7329749511860d-01, &
   0.9226316616509d-10, 0.1259158839583d+01, 0.1990721704425d+00, &
   0.8664946419555d-10, 0.3353244696934d+01, 0.5439178814476d+00, &
   0.7757230468978d-10, 0.1447677295196d+01, 0.9491756770005d+00, &
   0.7693168628139d-10, 0.1120509896721d+01, 0.1030928125552d+00, &
   0.5487897454612d-10, 0.4439380426795d+01, 0.8531963191132d+00 /
   DATA ((S1(I, J, 2), I = 1, 3), J = 31, NS1Y) /&
   0.5196118677218d-10, 0.3788856619137d+00, 0.2093666171530d+00, &
   0.5110853339935d-10, 0.1386879372016d+01, 0.2172315424036d+00, &
   0.5027804534813d-10, 0.1647881805466d+00, 0.2164800718209d+00, &
   0.4922485922674d-10, 0.1594315079862d+01, 0.2101180877357d+00, &
   0.6155599524400d-10, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.4447147832161d-10, 0.5480720918976d+01, 0.3235053470014d+00, &
   0.4144691276422d-10, 0.1931371033660d+01, 0.6373574839730d-01, &
   0.4099950625452d-10, 0.5229611294335d+01, 0.6470106940028d+00, &
   0.5060541682953d-10, 0.1731112486298d+01, 0.1422690933580d-01, &
   0.4293615946300d-10, 0.2714571038925d+01, 0.7358765972222d+00 /

   !  SS-to-Sun, T^2, Y
   DATA ((S2(I, J, 2), I = 1, 3), J = 1, NS2Y) /&
   0.1609114495091d-11, 0.2831096993481d+01, 0.2061856251104d+00, &
   0.1560330784946d-11, 0.5193058213906d+01, 0.2204125344462d+00, &
   0.1183535479202d-11, 0.5707003443890d+01, 0.5225775174439d+00, &
   0.1158183066182d-11, 0.1782400404928d+01, 0.5368044267797d+00, &
   0.1032868027407d-11, 0.4036925452011d+01, 0.2132990797783d+00, &
   0.6540142847741d-12, 0.4058241056717d+01, 0.4265981595566d+00, &
   0.7305236491596d-12, 0.6175401942957d+00, 0.5296909721118d+00, &
   -0.5580725052968d-12, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3946122651015d-12, 0.4108265279171d+00, 0.1059381944224d+01 /

   !  SS-to-Sun, T^0, Z
   DATA ((S0(I, J, 3), I = 1, 3), J = 1, 10) /&
   0.1181255122986d-03, 0.4607918989164d+00, 0.2132990797783d+00, &
   0.1127777651095d-03, 0.4169146331296d+00, 0.5296909721118d+00, &
   0.4777754401806d-04, 0.4582657007130d+01, 0.3813291813120d-01, &
   0.1129354285772d-04, 0.5758735142480d+01, 0.7478166569050d-01, &
   -0.1149543637123d-04, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3298730512306d-05, 0.5978801994625d+01, 0.4265981595566d+00, &
   0.2733376706079d-05, 0.7665413691040d+00, 0.1059381944224d+01, &
   0.9426389657270d-06, 0.3710201265838d+01, 0.2061856251104d+00, &
   0.8187517749552d-06, 0.3390675605802d+00, 0.2204125344462d+00, &
   0.4080447871819d-06, 0.4552296640088d+00, 0.5225775174439d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 11, 20) /&
   0.3169973017028d-06, 0.3445455899321d+01, 0.5368044267797d+00, &
   0.2438098615549d-06, 0.5664675150648d+01, 0.3664874755930d-01, &
   0.2601897517235d-06, 0.1931894095697d+01, 0.1495633313810d+00, &
   0.2314558080079d-06, 0.3666319115574d+00, 0.3961708870310d-01, &
   0.1962549548002d-06, 0.3167411699020d+01, 0.7626583626240d-01, &
   0.2180518287925d-06, 0.1544420746580d+01, 0.7113454667900d-02, &
   0.1451382442868d-06, 0.1583756740070d+01, 0.1102062672231d+00, &
   0.1358439007389d-06, 0.5239941758280d+01, 0.6398972393349d+00, &
   0.1050585898028d-06, 0.2266958352859d+01, 0.3163918923335d+00, &
   0.1050029870186d-06, 0.2711495250354d+01, 0.4194847048887d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 21, 30) /&
   0.9934920679800d-07, 0.1116208151396d+01, 0.1589072916335d+01, &
   0.1048395331560d-06, 0.3408619600206d+01, 0.1021328554739d+02, &
   0.8370147196668d-07, 0.3810459401087d+01, 0.2535050500000d-01, &
   0.7989856510998d-07, 0.3769910473647d+01, 0.7329749511860d-01, &
   0.5441221655233d-07, 0.2416994903374d+01, 0.1030928125552d+00, &
   0.4610812906784d-07, 0.5858503336994d+01, 0.4337116142245d+00, &
   0.3923022803444d-07, 0.3354170010125d+00, 0.1484170571900d-02, &
   0.2610725582128d-07, 0.5410600646324d+01, 0.6327837846670d+00, &
   0.2455279767721d-07, 0.6120216681403d+01, 0.1162474756779d+01, &
   0.2375530706525d-07, 0.6055443426143d+01, 0.1052268489556d+01 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 31, 40) /&
   0.1782967577553d-07, 0.3146108708004d+01, 0.8460828644453d+00, &
   0.1581687095238d-07, 0.6255496089819d+00, 0.3340612434717d+01, &
   0.1594657672461d-07, 0.3782604300261d+01, 0.1066495398892d+01, &
   0.1563448615040d-07, 0.1997775733196d+01, 0.2022531624851d+00, &
   0.1463624258525d-07, 0.1736316792088d+00, 0.3516457698740d-01, &
   0.1331585056673d-07, 0.4331941830747d+01, 0.9491756770005d+00, &
   0.1130634557637d-07, 0.6152017751825d+01, 0.2968341143800d-02, &
   0.1028949607145d-07, 0.2101792614637d+00, 0.2275259891141d+00, &
   0.1024074971618d-07, 0.4071833211074d+01, 0.5070101000000d-01, &
   0.8826956060303d-08, 0.4861633688145d+00, 0.2093666171530d+00 /
   DATA ((S0(I, J, 3), I = 1, 3), J = 41, NS0Z) /&
   0.8572230171541d-08, 0.5268190724302d+01, 0.4110125927500d-01, &
   0.7649332643544d-08, 0.5134543417106d+01, 0.2608790314060d+02, &
   0.8581673291033d-08, 0.2920218146681d+01, 0.1480791608091d+00, &
   0.8430589300938d-08, 0.3604576619108d+01, 0.2172315424036d+00, &
   0.7776165501012d-08, 0.3772942249792d+01, 0.6373574839730d-01, &
   0.8311070234408d-08, 0.6200412329888d+01, 0.3235053470014d+00, &
   0.6927365212582d-08, 0.4543353113437d+01, 0.8531963191132d+00, &
   0.6791574208598d-08, 0.2882188406238d+01, 0.7181332454670d-01, &
   0.5593100811839d-08, 0.1776646892780d+01, 0.7429900518901d+00, &
   0.7788777276590d-09, 0.1900569908215d+01, 0.5217580628120d+02, &
   0.4553381853021d-08, 0.3949617611240d+01, 0.7775000683430d-01 /

   !  SS-to-Sun, T^1, Z
   DATA ((S1(I, J, 3), I = 1, 3), J = 1, 10) /&
   0.5444220475678d-08, 0.1803825509310d+01, 0.2132990797783d+00, &
   0.3883412695596d-08, 0.4668616389392d+01, 0.5296909721118d+00, &
   0.1334341434551d-08, 0.0000000000000d+00, 0.0000000000000d+00, &
   0.3730001266883d-09, 0.5401405918943d+01, 0.2061856251104d+00, &
   0.2894929197956d-09, 0.4932415609852d+01, 0.2204125344462d+00, &
   0.2857950357701d-09, 0.3154625362131d+01, 0.7478166569050d-01, &
   0.2499226432292d-09, 0.3657486128988d+01, 0.4265981595566d+00, &
   0.1937705443593d-09, 0.5740434679002d+01, 0.1059381944224d+01, &
   0.1374894396320d-09, 0.1712857366891d+01, 0.5368044267797d+00, &
   0.1217248678408d-09, 0.2312090870932d+01, 0.5225775174439d+00 /
   DATA ((S1(I, J, 3), I = 1, 3), J = 11, NS1Z) /&
   0.7961052740870d-10, 0.5283368554163d+01, 0.3813291813120d-01, &
   0.4979225949689d-10, 0.4298290471860d+01, 0.4194847048887d+00, &
   0.4388552286597d-10, 0.6145515047406d+01, 0.7113454667900d-02, &
   0.2586835212560d-10, 0.3019448001809d+01, 0.6398972393349d+00 /

   !  SS-to-Sun, T^2, Z
   DATA ((S2(I, J, 3), I = 1, 3), J = 1, NS2Z) /&
   0.3749920358054d-12, 0.3230285558668d+01, 0.2132990797783d+00, &
   0.2735037220939d-12, 0.6154322683046d+01, 0.5296909721118d+00 /

   t = (xjd-xjd2000)/yearj
   t2= t * t
   do j = 1, 3 ! loop on x, y, z
      aa = 0d0
      bb = 0d0

      do k = 1, ns0(j)
         arg = s0(2, k, j) + s0(3, k, j) * t
         amp = s0(1, k, j)
         aa = aa + amp * cos(arg)
         bb = bb - s0(3, k, j) * amp * sin(arg)
      enddo

      do k = 1, ns1(j)
         arg = s1(2, k, j) + s1(3, k, j) * t
         cc = cos(arg)
         amp = s1(1, k, j)
         aa = aa + amp * t * cc
         bb = bb + (amp * cc - s1(3, k, j) * amp * t * sin(arg))
      enddo

      do k = 1, ns2(j)
         arg = s2(2, k, j) + s2(3, k, j) * t
         cc = cos(arg)
         amp = s2(1, k, j)
         aa = aa + amp * t2 * cc
         bb = bb + (2d0 * t * amp * cc - s2(3, k, j) * amp * t2 * sin(arg))
      enddo

      xyz(j) = aa ! position in au
      xyzp(j) = bb/yearj ! velocity in au/day
   enddo
   return
end subroutine baryc_hm
!********************************************************
subroutine baryc_h(xjd, xyz, xyzp)
   !********************************************************
   !
   !     Vector Barycenter to  Sun at the julian date xjd
   !     xyz(i)  in au,  xyzp in au/day
   !     ecliptic frame at J2000.
   !     Highest accuracy
   !
   !     Source ! from VSOP Series provided by P. Bretagnon June 2000
   !
   !     F. Mignard OCA/Cassiopee
   !
   !    *INPUT
   !     xjd     : julian day
   !
   !    *OUTPUT
   !     xyz     : Sun from SS barycenter in au in ecliptic frame J2000
   !     xyzp    : Velocity of the Sun in au/day
   !
   !     Precision : 1.5e-8  au      ~  2.5  km
   !                 3.0e-10 au/day  ~  0.5 mm/s
   !
   !********************************************************
   !

   implicit none

   integer, parameter :: ns0m = 213, ns1m = 50, ns2m = 9
   integer, parameter :: ns0x = 212, ns1x = 50, ns2x = 9
   integer, parameter :: ns0y = 213, ns1y = 50, ns2y = 9
   integer, parameter :: ns0z = 69, ns1z = 14, ns2z = 2
   
   

   real(kind = dp), intent(in) :: xjd
   real(kind = dp), dimension(3), intent(out) :: xyz, xyzp

   integer :: ns0(3) = (/ns0x, ns0y, ns0z/)
   integer :: ns1(3) = (/ns1x, ns1y, ns1z/)
   integer :: ns2(3) = (/ns2x, ns2y, ns2z/)

   real(kind = dp) :: s0(3, ns0m, 3), s1(3, ns1m, 3), s2(3, ns2m, 3)

   real(kind = dp) :: t, t2, aa, bb, arg, amp, cc
   integer :: i, j, k


  !  SSB-to-Sun, T^0, X
  DATA ((S0(I, J, 1), I = 1, 3), J = 1, 10) /&
  0.4956757536410d-02, 0.3741073751789d+01, 0.5296909721118d+00, &
  0.2718490072522d-02, 0.4016011511425d+01, 0.2132990797783d+00, &
  0.1546493974344d-02, 0.2170528330642d+01, 0.3813291813120d-01, &
  0.8366855276341d-03, 0.2339614075294d+01, 0.7478166569050d-01, &
  0.2936777942117d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.1201317439469d-03, 0.4090736353305d+01, 0.1059381944224d+01, &
  0.7578550887230d-04, 0.3241518088140d+01, 0.4265981595566d+00, &
  0.1941787367773d-04, 0.1012202064330d+01, 0.2061856251104d+00, &
  0.1889227765991d-04, 0.3892520416440d+01, 0.2204125344462d+00, &
  0.1937896968613d-04, 0.4797779441161d+01, 0.1495633313810d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 11, 20) /&
  0.1434506110873d-04, 0.3868960697933d+01, 0.5225775174439d+00, &
  0.1406659911580d-04, 0.4759766557397d+00, 0.5368044267797d+00, &
  0.1179022300202d-04, 0.7774961520598d+00, 0.7626583626240d-01, &
  0.8085864460959d-05, 0.3254654471465d+01, 0.3664874755930d-01, &
  0.7622752967615d-05, 0.4227633103489d+01, 0.3961708870310d-01, &
  0.6209171139066d-05, 0.2791828325711d+00, 0.7329749511860d-01, &
  0.4366435633970d-05, 0.4440454875925d+01, 0.1589072916335d+01, &
  0.3792124889348d-05, 0.5156393842356d+01, 0.7113454667900d-02, &
  0.3154548963402d-05, 0.6157005730093d+01, 0.4194847048887d+00, &
  0.3088359882942d-05, 0.2494567553163d+01, 0.6398972393349d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 21, 30) /&
  0.2788440902136d-05, 0.4934318747989d+01, 0.1102062672231d+00, &
  0.3039928456376d-05, 0.4895077702640d+01, 0.6283075850446d+01, &
  0.2272258457679d-05, 0.5278394064764d+01, 0.1030928125552d+00, &
  0.2162007057957d-05, 0.5802978019099d+01, 0.3163918923335d+00, &
  0.1767632855737d-05, 0.3415346595193d-01, 0.1021328554739d+02, &
  0.1349413459362d-05, 0.2001643230755d+01, 0.1484170571900d-02, &
  0.1170141900476d-05, 0.2424750491620d+01, 0.6327837846670d+00, &
  0.1054355266820d-05, 0.3123311487576d+01, 0.4337116142245d+00, &
  0.9800822461610d-06, 0.3026258088130d+01, 0.1052268489556d+01, &
  0.1091203749931d-05, 0.3157811670347d+01, 0.1162474756779d+01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 31, 40) /&
  0.6960236715913d-06, 0.8219570542313d+00, 0.1066495398892d+01, &
  0.5689257296909d-06, 0.1323052375236d+01, 0.9491756770005d+00, &
  0.6613172135802d-06, 0.2765348881598d+00, 0.8460828644453d+00, &
  0.6277702517571d-06, 0.5794064466382d+01, 0.1480791608091d+00, &
  0.6304884066699d-06, 0.7323555380787d+00, 0.2243449970715d+00, &
  0.4897850467382d-06, 0.3062464235399d+01, 0.3340612434717d+01, &
  0.3759148598786d-06, 0.4588290469664d+01, 0.3516457698740d-01, &
  0.3110520548195d-06, 0.1374299536572d+01, 0.6373574839730d-01, &
  0.3064708359780d-06, 0.4222267485047d+01, 0.1104591729320d-01, &
  0.2856347168241d-06, 0.3714202944973d+01, 0.1510475019529d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 41, 50) /&
  0.2840945514288d-06, 0.2847972875882d+01, 0.4110125927500d-01, &
  0.2378951599405d-06, 0.3762072563388d+01, 0.2275259891141d+00, &
  0.2714229481417d-06, 0.1036049980031d+01, 0.2535050500000d-01, &
  0.2323551717307d-06, 0.4682388599076d+00, 0.8582758298370d-01, &
  0.1881790512219d-06, 0.4790565425418d+01, 0.2118763888447d+01, &
  0.2261353968371d-06, 0.1669144912212d+01, 0.7181332454670d-01, &
  0.2214546389848d-06, 0.3937717281614d+01, 0.2968341143800d-02, &
  0.2184915594933d-06, 0.1129169845099d+00, 0.7775000683430d-01, &
  0.2000164937936d-06, 0.4030009638488d+01, 0.2093666171530d+00, &
  0.1966105136719d-06, 0.8745955786834d+00, 0.2172315424036d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 51, 60) /&
  0.1904742332624d-06, 0.5919743598964d+01, 0.2022531624851d+00, &
  0.1657399705031d-06, 0.2549141484884d+01, 0.7358765972222d+00, &
  0.1574070533987d-06, 0.5277533020230d+01, 0.7429900518901d+00, &
  0.1832261651039d-06, 0.3064688127777d+01, 0.3235053470014d+00, &
  0.1733615346569d-06, 0.3011432799094d+01, 0.1385174140878d+00, &
  0.1549124014496d-06, 0.4005569132359d+01, 0.5154640627760d+00, &
  0.1637044713838d-06, 0.1831375966632d+01, 0.8531963191132d+00, &
  0.1123420082383d-06, 0.1180270407578d+01, 0.1990721704425d+00, &
  0.1083754165740d-06, 0.3414101320863d+00, 0.5439178814476d+00, &
  0.1156638012655d-06, 0.6130479452594d+00, 0.5257585094865d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 61, 70) /&
  0.1142548785134d-06, 0.3724761948846d+01, 0.5336234347371d+00, &
  0.7921463895965d-07, 0.2435425589361d+01, 0.1478866649112d+01, &
  0.7428600285231d-07, 0.3542144398753d+01, 0.2164800718209d+00, &
  0.8323211246747d-07, 0.3525058072354d+01, 0.1692165728891d+01, &
  0.7257595116312d-07, 0.1364299431982d+01, 0.2101180877357d+00, &
  0.7111185833236d-07, 0.2460478875808d+01, 0.4155522422634d+00, &
  0.6868090383716d-07, 0.4397327670704d+01, 0.1173197218910d+00, &
  0.7226419974175d-07, 0.4042647308905d+01, 0.1265567569334d+01, &
  0.6955642383177d-07, 0.2865047906085d+01, 0.9562891316684d+00, &
  0.7492139296331d-07, 0.5014278994215d+01, 0.1422690933580d-01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 71, 80) /&
  0.6598363128857d-07, 0.2376730020492d+01, 0.6470106940028d+00, &
  0.7381147293385d-07, 0.3272990384244d+01, 0.1581959461667d+01, &
  0.6402909624032d-07, 0.5302290955138d+01, 0.9597935788730d-01, &
  0.6237454263857d-07, 0.5444144425332d+01, 0.7084920306520d-01, &
  0.5241198544016d-07, 0.4215359579205d+01, 0.5265099800692d+00, &
  0.5144463853918d-07, 0.1218916689916d+00, 0.5328719641544d+00, &
  0.5868164772299d-07, 0.2369402002213d+01, 0.7871412831580d-01, &
  0.6233195669151d-07, 0.1254922242403d+01, 0.2608790314060d+02, &
  0.6068463791422d-07, 0.5679713760431d+01, 0.1114304132498d+00, &
  0.4359361135065d-07, 0.6097219641646d+00, 0.1375773836557d+01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 81, 90) /&
  0.4686510366826d-07, 0.4786231041431d+01, 0.1143987543936d+00, &
  0.3758977287225d-07, 0.1167368068139d+01, 0.1596186371003d+01, &
  0.4282051974778d-07, 0.1519471064319d+01, 0.2770348281756d+00, &
  0.5153765386113d-07, 0.1860532322984d+01, 0.2228608264996d+00, &
  0.4575129387188d-07, 0.7632857887158d+00, 0.1465949902372d+00, &
  0.3326844933286d-07, 0.1298219485285d+01, 0.5070101000000d-01, &
  0.3748617450984d-07, 0.1046510321062d+01, 0.4903339079539d+00, &
  0.2816756661499d-07, 0.3434522346190d+01, 0.2991266627620d+00, &
  0.3412750405039d-07, 0.2523766270318d+01, 0.3518164938661d+00, &
  0.2655796761776d-07, 0.2904422260194d+01, 0.6256703299991d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 91, 100) /&
  0.2963597929458d-07, 0.5923900431149d+00, 0.1099462426779d+00, &
  0.2539523734781d-07, 0.4851947722567d+01, 0.1256615170089d+02, &
  0.2283087914139d-07, 0.3400498595496d+01, 0.6681224869435d+01, &
  0.2321309799331d-07, 0.5789099148673d+01, 0.3368040641550d-01, &
  0.2549657649750d-07, 0.3991856479792d-01, 0.1169588211447d+01, &
  0.2290462303977d-07, 0.2788567577052d+01, 0.1045155034888d+01, &
  0.1945398522914d-07, 0.3290896998176d+01, 0.1155361302111d+01, &
  0.1849171512638d-07, 0.2698060129367d+01, 0.4452511715700d-02, &
  0.1647199834254d-07, 0.3016735644085d+01, 0.4408250688924d+00, &
  0.1529530765273d-07, 0.5573043116178d+01, 0.6521991896920d-01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 101, 110) /&
  0.1433199339978d-07, 0.1481192356147d+01, 0.9420622223326d+00, &
  0.1729134193602d-07, 0.1422817538933d+01, 0.2108507877249d+00, &
  0.1716463931346d-07, 0.3469468901855d+01, 0.2157473718317d+00, &
  0.1391206061378d-07, 0.6122436220547d+01, 0.4123712502208d+00, &
  0.1404746661924d-07, 0.1647765641936d+01, 0.4258542984690d-01, &
  0.1410452399455d-07, 0.5989729161964d+01, 0.2258291676434d+00, &
  0.1089828772168d-07, 0.2833705509371d+01, 0.4226656969313d+00, &
  0.1047374564948d-07, 0.5090690007331d+00, 0.3092784376656d+00, &
  0.1358279126532d-07, 0.5128990262836d+01, 0.7923417740620d-01, &
  0.1020456476148d-07, 0.9632772880808d+00, 0.1456308687557d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 111, 120) /&
  0.1033428735328d-07, 0.3223779318418d+01, 0.1795258541446d+01, &
  0.1412435841540d-07, 0.2410271572721d+01, 0.1525316725248d+00, &
  0.9722759371574d-08, 0.2333531395690d+01, 0.8434341241180d-01, &
  0.9657334084704d-08, 0.6199270974168d+01, 0.1272681024002d+01, &
  0.1083641148690d-07, 0.2864222292929d+01, 0.7032915397480d-01, &
  0.1067318403838d-07, 0.5833458866568d+00, 0.2123349582968d+00, &
  0.1062366201976d-07, 0.4307753989494d+01, 0.2142632012598d+00, &
  0.1236364149266d-07, 0.2873917870593d+01, 0.1847279083684d+00, &
  0.1092759489593d-07, 0.2959887266733d+01, 0.1370332435159d+00, &
  0.8912069362899d-08, 0.5141213702562d+01, 0.2648454860559d+01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 121, 130) /&
  0.9656467707970d-08, 0.4532182462323d+01, 0.4376440768498d+00, &
  0.8098386150135d-08, 0.2268906338379d+01, 0.2880807454688d+00, &
  0.7857714675000d-08, 0.4055544260745d+01, 0.2037373330570d+00, &
  0.7288455940646d-08, 0.5357901655142d+01, 0.1129145838217d+00, &
  0.9450595950552d-08, 0.4264926963939d+01, 0.5272426800584d+00, &
  0.9381718247537d-08, 0.7489366976576d-01, 0.5321392641652d+00, &
  0.7079052646038d-08, 0.1923311052874d+01, 0.6288513220417d+00, &
  0.9259004415344d-08, 0.2970256853438d+01, 0.1606092486742d+00, &
  0.8259801499742d-08, 0.3327056314697d+01, 0.8389694097774d+00, &
  0.6476334355779d-08, 0.2954925505727d+01, 0.2008557621224d+01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 131, 140) /&
  0.5984021492007d-08, 0.9138753105829d+00, 0.2042657109477d+02, &
  0.5989546863181d-08, 0.3244464082031d+01, 0.2111650433779d+01, &
  0.6233108606023d-08, 0.4995232638403d+00, 0.4305306221819d+00, &
  0.6877299149965d-08, 0.2834987233449d+01, 0.9561746721300d-02, &
  0.8311234227190d-08, 0.2202951835758d+01, 0.3801276407308d+00, &
  0.6599472832414d-08, 0.4478581462618d+01, 0.1063314406849d+01, &
  0.6160491096549d-08, 0.5145858696411d+01, 0.1368660381889d+01, &
  0.6164772043891d-08, 0.3762976697911d+00, 0.4234171675140d+00, &
  0.6363248684450d-08, 0.3162246718685d+01, 0.1253008786510d-01, &
  0.6448587520999d-08, 0.3442693302119d+01, 0.5287268506303d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 141, 150) /&
  0.6431662283977d-08, 0.8977549136606d+00, 0.5306550935933d+00, &
  0.6351223158474d-08, 0.4306447410369d+01, 0.5217580628120d+02, &
  0.5476721393451d-08, 0.3888529177855d+01, 0.2221856701002d+01, &
  0.5341772572619d-08, 0.2655560662512d+01, 0.7466759693650d-01, &
  0.5337055758302d-08, 0.5164990735946d+01, 0.7489573444450d-01, &
  0.5373120816787d-08, 0.6041214553456d+01, 0.1274714967946d+00, &
  0.5392351705426d-08, 0.9177763485932d+00, 0.1055449481598d+01, &
  0.6688495850205d-08, 0.3089608126937d+01, 0.2213766559277d+00, &
  0.5072003660362d-08, 0.4311316541553d+01, 0.2132517061319d+00, &
  0.5070726650455d-08, 0.5790675464444d+00, 0.2133464534247d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 151, 160) /&
  0.5658012950032d-08, 0.2703945510675d+01, 0.7287631425543d+00, &
  0.4835509924854d-08, 0.2975422976065d+01, 0.7160067364790d-01, &
  0.6479821978012d-08, 0.1324168733114d+01, 0.2209183458640d-01, &
  0.6230636494980d-08, 0.2860103632836d+01, 0.3306188016693d+00, &
  0.4649239516213d-08, 0.4832259763403d+01, 0.7796265773310d-01, &
  0.6487325792700d-08, 0.2726165825042d+01, 0.3884652414254d+00, &
  0.4682823682770d-08, 0.6966602455408d+00, 0.1073608853559d+01, &
  0.5704230804976d-08, 0.5669634104606d+01, 0.8731175355560d-01, &
  0.6125413585489d-08, 0.1513386538915d+01, 0.7605151500000d-01, &
  0.6035825038187d-08, 0.1983509168227d+01, 0.9846002785331d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 161, 170) /&
  0.4331123462303d-08, 0.2782892992807d+01, 0.4297791515992d+00, &
  0.4681107685143d-08, 0.5337232886836d+01, 0.2127790306879d+00, &
  0.4669105829655d-08, 0.5837133792160d+01, 0.2138191288687d+00, &
  0.5138823602365d-08, 0.3080560200507d+01, 0.7233337363710d-01, &
  0.4615856664534d-08, 0.1661747897471d+01, 0.8603097737811d+00, &
  0.4496916702197d-08, 0.2112508027068d+01, 0.7381754420900d-01, &
  0.4278479042945d-08, 0.5716528462627d+01, 0.7574578717200d-01, &
  0.3840525503932d-08, 0.6424172726492d+00, 0.3407705765729d+00, &
  0.4866636509685d-08, 0.4919244697715d+01, 0.7722995774390d-01, &
  0.3526100639296d-08, 0.2550821052734d+01, 0.6225157782540d-01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 171, 180) /&
  0.3939558488075d-08, 0.3939331491710d+01, 0.5268983110410d-01, &
  0.4041268772576d-08, 0.2275337571218d+01, 0.3503323232942d+00, &
  0.3948761842853d-08, 0.1999324200790d+01, 0.1451108196653d+00, &
  0.3258394550029d-08, 0.9121001378200d+00, 0.5296435984654d+00, &
  0.3257897048761d-08, 0.3428428660869d+01, 0.5297383457582d+00, &
  0.3842559031298d-08, 0.6132927720035d+01, 0.9098186128426d+00, &
  0.3109920095448d-08, 0.7693650193003d+00, 0.3932462625300d-02, &
  0.3132237775119d-08, 0.3621293854908d+01, 0.2346394437820d+00, &
  0.3942189421510d-08, 0.4841863659733d+01, 0.3180992042600d-02, &
  0.3796972285340d-08, 0.1814174994268d+01, 0.1862120789403d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 181, 190) /&
  0.3995640233688d-08, 0.1386990406091d+01, 0.4549093064213d+00, &
  0.2875013727414d-08, 0.9178318587177d+00, 0.1905464808669d+01, &
  0.3073719932844d-08, 0.2688923811835d+01, 0.3628624111593d+00, &
  0.2731016580075d-08, 0.1188259127584d+01, 0.2131850110243d+00, &
  0.2729549896546d-08, 0.3702160634273d+01, 0.2134131485323d+00, &
  0.3339372892449d-08, 0.7199163960331d+00, 0.2007689919132d+00, &
  0.2898833764204d-08, 0.1916709364999d+01, 0.5291709230214d+00, &
  0.2894536549362d-08, 0.2424043195547d+01, 0.5302110212022d+00, &
  0.3096872473843d-08, 0.4445894977497d+01, 0.2976424921901d+00, &
  0.2635672326810d-08, 0.3814366984117d+01, 0.1485980103780d+01 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 191, 200) /&
  0.3649302697001d-08, 0.2924200596084d+01, 0.6044726378023d+00, &
  0.3127954585895d-08, 0.1842251648327d+01, 0.1084620721060d+00, &
  0.2616040173947d-08, 0.4155841921984d+01, 0.1258454114666d+01, &
  0.2597395859860d-08, 0.1158045978874d+00, 0.2103781122809d+00, &
  0.2593286172210d-08, 0.4771850408691d+01, 0.2162200472757d+00, &
  0.2481823585747d-08, 0.4608842558889d+00, 0.1062562936266d+01, &
  0.2742219550725d-08, 0.1538781127028d+01, 0.5651155736444d+00, &
  0.3199558469610d-08, 0.3226647822878d+00, 0.7036329877322d+00, &
  0.2666088542957d-08, 0.1967991731219d+00, 0.1400015846597d+00, &
  0.2397067430580d-08, 0.3707036669873d+01, 0.2125476091956d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 201, 210) /&
  0.2376570772738d-08, 0.1182086628042d+01, 0.2140505503610d+00, &
  0.2547228007887d-08, 0.4906256820629d+01, 0.1534957940063d+00, &
  0.2265575594114d-08, 0.3414949866857d+01, 0.2235935264888d+00, &
  0.2464381430585d-08, 0.4599122275378d+01, 0.2091065926078d+00, &
  0.2433408527044d-08, 0.2830751145445d+00, 0.2174915669488d+00, &
  0.2443605509076d-08, 0.4212046432538d+01, 0.1739420156204d+00, &
  0.2319779262465d-08, 0.9881978408630d+00, 0.7530171478090d-01, &
  0.2284622835465d-08, 0.5565347331588d+00, 0.7426161660010d-01, &
  0.2467268750783d-08, 0.5655708150766d+00, 0.2526561439362d+00, &
  0.2808513492782d-08, 0.1418405053408d+01, 0.5636314030725d+00 /
  DATA ((S0(I, J, 1), I = 1, 3), J = 211, NS0X) /&
  0.2329528932532d-08, 0.4069557545675d+01, 0.1056200952181d+01, &
  0.9698639532817d-09, 0.1074134313634d+01, 0.7826370942180d+02 /

  !  SSB-to-Sun, T^1, X
  DATA ((S1(I, J, 1), I = 1, 3), J = 1, 10) /&
  -0.1296310361520d-07, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.8975769009438d-08, 0.1128891609250d+01, 0.4265981595566d+00, &
  0.7771113441307d-08, 0.2706039877077d+01, 0.2061856251104d+00, &
  0.7538303866642d-08, 0.2191281289498d+01, 0.2204125344462d+00, &
  0.6061384579336d-08, 0.3248167319958d+01, 0.1059381944224d+01, &
  0.5726994235594d-08, 0.5569981398610d+01, 0.5225775174439d+00, &
  0.5616492836424d-08, 0.5057386614909d+01, 0.5368044267797d+00, &
  0.1010881584769d-08, 0.3473577116095d+01, 0.7113454667900d-02, &
  0.7259606157626d-09, 0.3651858593665d+00, 0.6398972393349d+00, &
  0.8755095026935d-09, 0.1662835408338d+01, 0.4194847048887d+00 /
  DATA ((S1(I, J, 1), I = 1, 3), J = 11, 20) /&
  0.5370491182812d-09, 0.1327673878077d+01, 0.4337116142245d+00, &
  0.5743773887665d-09, 0.4250200846687d+01, 0.2132990797783d+00, &
  0.4408103140300d-09, 0.3598752574277d+01, 0.1589072916335d+01, &
  0.3101892374445d-09, 0.4887822983319d+01, 0.1052268489556d+01, &
  0.3209453713578d-09, 0.9702272295114d+00, 0.5296909721118d+00, &
  0.3017228286064d-09, 0.5484462275949d+01, 0.1066495398892d+01, &
  0.3200700038601d-09, 0.2846613338643d+01, 0.1495633313810d+00, &
  0.2137637279911d-09, 0.5692163292729d+00, 0.3163918923335d+00, &
  0.1899686386727d-09, 0.2061077157189d+01, 0.2275259891141d+00, &
  0.1401994545308d-09, 0.4177771136967d+01, 0.1102062672231d+00 /
  DATA ((S1(I, J, 1), I = 1, 3), J = 21, 30) /&
  0.1578057810499d-09, 0.5782460597335d+01, 0.7626583626240d-01, &
  0.1237713253351d-09, 0.5705900866881d+01, 0.5154640627760d+00, &
  0.1313076837395d-09, 0.5163438179576d+01, 0.3664874755930d-01, &
  0.1184963304860d-09, 0.3054804427242d+01, 0.6327837846670d+00, &
  0.1238130878565d-09, 0.2317292575962d+01, 0.3961708870310d-01, &
  0.1015959527736d-09, 0.2194643645526d+01, 0.7329749511860d-01, &
  0.9017954423714d-10, 0.2868603545435d+01, 0.1990721704425d+00, &
  0.8668024955603d-10, 0.4923849675082d+01, 0.5439178814476d+00, &
  0.7756083930103d-10, 0.3014334135200d+01, 0.9491756770005d+00, &
  0.7536503401741d-10, 0.2704886279769d+01, 0.1030928125552d+00 /
  DATA ((S1(I, J, 1), I = 1, 3), J = 31, 40) /&
  0.5483308679332d-10, 0.6010983673799d+01, 0.8531963191132d+00, &
  0.5184339620428d-10, 0.1952704573291d+01, 0.2093666171530d+00, &
  0.5108658712030d-10, 0.2958575786649d+01, 0.2172315424036d+00, &
  0.5019424524650d-10, 0.1736317621318d+01, 0.2164800718209d+00, &
  0.4909312625978d-10, 0.3167216416257d+01, 0.2101180877357d+00, &
  0.4456638901107d-10, 0.7697579923471d+00, 0.3235053470014d+00, &
  0.4227030350925d-10, 0.3490910137928d+01, 0.6373574839730d-01, &
  0.4095456040093d-10, 0.5178888984491d+00, 0.6470106940028d+00, &
  0.4990537041422d-10, 0.3323887668974d+01, 0.1422690933580d-01, &
  0.4321170010845d-10, 0.4288484987118d+01, 0.7358765972222d+00 /
  DATA ((S1(I, J, 1), I = 1, 3), J = 41, NS1X) /&
  0.3544072091802d-10, 0.6021051579251d+01, 0.5265099800692d+00, &
  0.3480198638687d-10, 0.4600027054714d+01, 0.5328719641544d+00, &
  0.3440287244435d-10, 0.4349525970742d+01, 0.8582758298370d-01, &
  0.3330628322713d-10, 0.2347391505082d+01, 0.1104591729320d-01, &
  0.2973060707184d-10, 0.4789409286400d+01, 0.5257585094865d+00, &
  0.2932606766089d-10, 0.5831693799927d+01, 0.5336234347371d+00, &
  0.2876972310953d-10, 0.2692638514771d+01, 0.1173197218910d+00, &
  0.2827488278556d-10, 0.2056052487960d+01, 0.2022531624851d+00, &
  0.2515028239756d-10, 0.7411863262449d+00, 0.9597935788730d-01, &
  0.2853033744415d-10, 0.3948481024894d+01, 0.2118763888447d+01 /

  !  SSB-to-Sun, T^2, X
  DATA ((S2(I, J, 1), I = 1, 3), J = 1, NS2X) /&
  0.1603551636587d-11, 0.4404109410481d+01, 0.2061856251104d+00, &
  0.1556935889384d-11, 0.4818040873603d+00, 0.2204125344462d+00, &
  0.1182594414915d-11, 0.9935762734472d+00, 0.5225775174439d+00, &
  0.1158794583180d-11, 0.3353180966450d+01, 0.5368044267797d+00, &
  0.9597358943932d-12, 0.5567045358298d+01, 0.2132990797783d+00, &
  0.6511516579605d-12, 0.5630872420788d+01, 0.4265981595566d+00, &
  0.7419792747688d-12, 0.2156188581957d+01, 0.5296909721118d+00, &
  0.3951972655848d-12, 0.1981022541805d+01, 0.1059381944224d+01, &
  0.4478223877045d-12, 0.0000000000000d+00, 0.0000000000000d+00 /

  !  SSB-to-Sun, T^0, Y
  DATA ((S0(I, J, 2), I = 1, 3), J = 1, 10) /&
  0.4955392320126d-02, 0.2170467313679d+01, 0.5296909721118d+00, &
  0.2722325167392d-02, 0.2444433682196d+01, 0.2132990797783d+00, &
  0.1546579925346d-02, 0.5992779281546d+00, 0.3813291813120d-01, &
  0.8363140252966d-03, 0.7687356310801d+00, 0.7478166569050d-01, &
  0.3385792683603d-03, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.1201192221613d-03, 0.2520035601514d+01, 0.1059381944224d+01, &
  0.7587125720554d-04, 0.1669954006449d+01, 0.4265981595566d+00, &
  0.1964155361250d-04, 0.5707743963343d+01, 0.2061856251104d+00, &
  0.1891900364909d-04, 0.2320960679937d+01, 0.2204125344462d+00, &
  0.1937373433356d-04, 0.3226940689555d+01, 0.1495633313810d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 11, 20) /&
  0.1437139941351d-04, 0.2301626908096d+01, 0.5225775174439d+00, &
  0.1406267683099d-04, 0.5188579265542d+01, 0.5368044267797d+00, &
  0.1178703080346d-04, 0.5489483248476d+01, 0.7626583626240d-01, &
  0.8079835186041d-05, 0.1683751835264d+01, 0.3664874755930d-01, &
  0.7623253594652d-05, 0.2656400462961d+01, 0.3961708870310d-01, &
  0.6248667483971d-05, 0.4992775362055d+01, 0.7329749511860d-01, &
  0.4366353695038d-05, 0.2869706279678d+01, 0.1589072916335d+01, &
  0.3829101568895d-05, 0.3572131359950d+01, 0.7113454667900d-02, &
  0.3175733773908d-05, 0.4535372530045d+01, 0.4194847048887d+00, &
  0.3092437902159d-05, 0.9230153317909d+00, 0.6398972393349d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 21, 30) /&
  0.2874168812154d-05, 0.3363143761101d+01, 0.1102062672231d+00, &
  0.3040119321826d-05, 0.3324250895675d+01, 0.6283075850446d+01, &
  0.2699723308006d-05, 0.2917882441928d+00, 0.1030928125552d+00, &
  0.2134832683534d-05, 0.4220997202487d+01, 0.3163918923335d+00, &
  0.1770412139433d-05, 0.4747318496462d+01, 0.1021328554739d+02, &
  0.1377264209373d-05, 0.4305058462401d+00, 0.1484170571900d-02, &
  0.1127814538960d-05, 0.8538177240740d+00, 0.6327837846670d+00, &
  0.1055608090130d-05, 0.1551800742580d+01, 0.4337116142245d+00, &
  0.9802673861420d-06, 0.1459646735377d+01, 0.1052268489556d+01, &
  0.1090329461951d-05, 0.1587351228711d+01, 0.1162474756779d+01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 31, 40) /&
  0.6959590025090d-06, 0.5534442628766d+01, 0.1066495398892d+01, &
  0.5664914529542d-06, 0.6030673003297d+01, 0.9491756770005d+00, &
  0.6607787763599d-06, 0.4989507233927d+01, 0.8460828644453d+00, &
  0.6269725742838d-06, 0.4222951804572d+01, 0.1480791608091d+00, &
  0.6301889697863d-06, 0.5444316669126d+01, 0.2243449970715d+00, &
  0.4891042662861d-06, 0.1490552839784d+01, 0.3340612434717d+01, &
  0.3457083123290d-06, 0.3030475486049d+01, 0.3516457698740d-01, &
  0.3032559967314d-06, 0.2652038793632d+01, 0.1104591729320d-01, &
  0.2841133988903d-06, 0.1276744786829d+01, 0.4110125927500d-01, &
  0.2855564444432d-06, 0.2143368674733d+01, 0.1510475019529d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 41, 50) /&
  0.2765157135038d-06, 0.5444186109077d+01, 0.6373574839730d-01, &
  0.2382312465034d-06, 0.2190521137593d+01, 0.2275259891141d+00, &
  0.2808060365077d-06, 0.5735195064841d+01, 0.2535050500000d-01, &
  0.2332175234405d-06, 0.9481985524859d-01, 0.7181332454670d-01, &
  0.2322488199659d-06, 0.5180499361533d+01, 0.8582758298370d-01, &
  0.1881850258423d-06, 0.3219788273885d+01, 0.2118763888447d+01, &
  0.2196111392808d-06, 0.2366941159761d+01, 0.2968341143800d-02, &
  0.2183810335519d-06, 0.4825445110915d+01, 0.7775000683430d-01, &
  0.2002733093326d-06, 0.2457148995307d+01, 0.2093666171530d+00, &
  0.1967111767229d-06, 0.5586291545459d+01, 0.2172315424036d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 51, 60) /&
  0.1568473250543d-06, 0.3708003123320d+01, 0.7429900518901d+00, &
  0.1852528314300d-06, 0.4310638151560d+01, 0.2022531624851d+00, &
  0.1832111226447d-06, 0.1494665322656d+01, 0.3235053470014d+00, &
  0.1746805502310d-06, 0.1451378500784d+01, 0.1385174140878d+00, &
  0.1555730966650d-06, 0.1068040418198d+01, 0.7358765972222d+00, &
  0.1554883462559d-06, 0.2442579035461d+01, 0.5154640627760d+00, &
  0.1638380568746d-06, 0.2597913420625d+00, 0.8531963191132d+00, &
  0.1159938593640d-06, 0.5834512021280d+01, 0.1990721704425d+00, &
  0.1083427965695d-06, 0.5054033177950d+01, 0.5439178814476d+00, &
  0.1156480369431d-06, 0.5325677432457d+01, 0.5257585094865d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 61, 70) /&
  0.1141308860095d-06, 0.2153403923857d+01, 0.5336234347371d+00, &
  0.7913146470946d-07, 0.8642846847027d+00, 0.1478866649112d+01, &
  0.7439752463733d-07, 0.1970628496213d+01, 0.2164800718209d+00, &
  0.7280277104079d-07, 0.6073307250609d+01, 0.2101180877357d+00, &
  0.8319567719136d-07, 0.1954371928334d+01, 0.1692165728891d+01, &
  0.7137705549290d-07, 0.8904989440909d+00, 0.4155522422634d+00, &
  0.6900825396225d-07, 0.2825717714977d+01, 0.1173197218910d+00, &
  0.7245757216635d-07, 0.2481677513331d+01, 0.1265567569334d+01, &
  0.6961165696255d-07, 0.1292955312978d+01, 0.9562891316684d+00, &
  0.7571804456890d-07, 0.3427517575069d+01, 0.1422690933580d-01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 71, 80) /&
  0.6605425721904d-07, 0.8052192701492d+00, 0.6470106940028d+00, &
  0.7375477357248d-07, 0.1705076390088d+01, 0.1581959461667d+01, &
  0.7041664951470d-07, 0.4848356967891d+00, 0.9597935788730d-01, &
  0.6322199535763d-07, 0.3878069473909d+01, 0.7084920306520d-01, &
  0.5244380279191d-07, 0.2645560544125d+01, 0.5265099800692d+00, &
  0.5143125704988d-07, 0.4834486101370d+01, 0.5328719641544d+00, &
  0.5871866319373d-07, 0.7981472548900d+00, 0.7871412831580d-01, &
  0.6300822573871d-07, 0.5979398788281d+01, 0.2608790314060d+02, &
  0.6062154271548d-07, 0.4108655402756d+01, 0.1114304132498d+00, &
  0.4361912339976d-07, 0.5322624319280d+01, 0.1375773836557d+01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 81, 90) /&
  0.4417005920067d-07, 0.6240817359284d+01, 0.2770348281756d+00, &
  0.4686806749936d-07, 0.3214977301156d+01, 0.1143987543936d+00, &
  0.3758892132305d-07, 0.5879809634765d+01, 0.1596186371003d+01, &
  0.5151351332319d-07, 0.2893377688007d+00, 0.2228608264996d+00, &
  0.4554683578572d-07, 0.5475427144122d+01, 0.1465949902372d+00, &
  0.3442381385338d-07, 0.5992034796640d+01, 0.5070101000000d-01, &
  0.2831093954933d-07, 0.5367350273914d+01, 0.3092784376656d+00, &
  0.3756267090084d-07, 0.5758171285420d+01, 0.4903339079539d+00, &
  0.2816374679892d-07, 0.1863718700923d+01, 0.2991266627620d+00, &
  0.3419307025569d-07, 0.9524347534130d+00, 0.3518164938661d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 91, 100) /&
  0.2904250494239d-07, 0.5304471615602d+01, 0.1099462426779d+00, &
  0.2471734511206d-07, 0.1297069793530d+01, 0.6256703299991d+00, &
  0.2539620831872d-07, 0.3281126083375d+01, 0.1256615170089d+02, &
  0.2281017868007d-07, 0.1829122133165d+01, 0.6681224869435d+01, &
  0.2275319473335d-07, 0.5797198160181d+01, 0.3932462625300d-02, &
  0.2547755368442d-07, 0.4752697708330d+01, 0.1169588211447d+01, &
  0.2285979669317d-07, 0.1223205292886d+01, 0.1045155034888d+01, &
  0.1913386560994d-07, 0.1757532993389d+01, 0.1155361302111d+01, &
  0.1809020525147d-07, 0.4246116108791d+01, 0.3368040641550d-01, &
  0.1649213300201d-07, 0.1445162890627d+01, 0.4408250688924d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 101, 110) /&
  0.1834972793932d-07, 0.1126917567225d+01, 0.4452511715700d-02, &
  0.1439550648138d-07, 0.6160756834764d+01, 0.9420622223326d+00, &
  0.1487645457041d-07, 0.4358761931792d+01, 0.4123712502208d+00, &
  0.1731729516660d-07, 0.6134456753344d+01, 0.2108507877249d+00, &
  0.1717747163567d-07, 0.1898186084455d+01, 0.2157473718317d+00, &
  0.1418190430374d-07, 0.4180286741266d+01, 0.6521991896920d-01, &
  0.1404844134873d-07, 0.7654053565412d-01, 0.4258542984690d-01, &
  0.1409842846538d-07, 0.4418612420312d+01, 0.2258291676434d+00, &
  0.1090948346291d-07, 0.1260615686131d+01, 0.4226656969313d+00, &
  0.1357577323612d-07, 0.3558248818690d+01, 0.7923417740620d-01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 111, 120) /&
  0.1018154061960d-07, 0.5676087241256d+01, 0.1456308687557d+00, &
  0.1412073972109d-07, 0.8394392632422d+00, 0.1525316725248d+00, &
  0.1030938326496d-07, 0.1653593274064d+01, 0.1795258541446d+01, &
  0.1180081567104d-07, 0.1285802592036d+01, 0.7032915397480d-01, &
  0.9708510575650d-08, 0.7631889488106d+00, 0.8434341241180d-01, &
  0.9637689663447d-08, 0.4630642649176d+01, 0.1272681024002d+01, &
  0.1068910429389d-07, 0.5294934032165d+01, 0.2123349582968d+00, &
  0.1063716179336d-07, 0.2736266800832d+01, 0.2142632012598d+00, &
  0.1234858713814d-07, 0.1302891146570d+01, 0.1847279083684d+00, &
  0.8912631189738d-08, 0.3570415993621d+01, 0.2648454860559d+01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 121, 130) /&
  0.1036378285534d-07, 0.4236693440949d+01, 0.1370332435159d+00, &
  0.9667798501561d-08, 0.2960768892398d+01, 0.4376440768498d+00, &
  0.8108314201902d-08, 0.6987781646841d+00, 0.2880807454688d+00, &
  0.7648364324628d-08, 0.2499017863863d+01, 0.2037373330570d+00, &
  0.7286136828406d-08, 0.3787426951665d+01, 0.1129145838217d+00, &
  0.9448237743913d-08, 0.2694354332983d+01, 0.5272426800584d+00, &
  0.9374276106428d-08, 0.4787121277064d+01, 0.5321392641652d+00, &
  0.7100226287462d-08, 0.3530238792101d+00, 0.6288513220417d+00, &
  0.9253056659571d-08, 0.1399478925664d+01, 0.1606092486742d+00, &
  0.6636432145504d-08, 0.3479575438447d+01, 0.1368660381889d+01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 131, 140) /&
  0.6469975312932d-08, 0.1383669964800d+01, 0.2008557621224d+01, &
  0.7335849729765d-08, 0.1243698166898d+01, 0.9561746721300d-02, &
  0.8743421205855d-08, 0.3776164289301d+01, 0.3801276407308d+00, &
  0.5993635744494d-08, 0.5627122113596d+01, 0.2042657109477d+02, &
  0.5981008479693d-08, 0.1674336636752d+01, 0.2111650433779d+01, &
  0.6188535145838d-08, 0.5214925208672d+01, 0.4305306221819d+00, &
  0.6596074017566d-08, 0.2907653268124d+01, 0.1063314406849d+01, &
  0.6630815126226d-08, 0.2127643669658d+01, 0.8389694097774d+00, &
  0.6156772830040d-08, 0.5082160803295d+01, 0.4234171675140d+00, &
  0.6446960563014d-08, 0.1872100916905d+01, 0.5287268506303d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 141, 150) /&
  0.6429324424668d-08, 0.5610276103577d+01, 0.5306550935933d+00, &
  0.6302232396465d-08, 0.1592152049607d+01, 0.1253008786510d-01, &
  0.6399244436159d-08, 0.2746214421532d+01, 0.5217580628120d+02, &
  0.5474965172558d-08, 0.2317666374383d+01, 0.2221856701002d+01, &
  0.5339293190692d-08, 0.1084724961156d+01, 0.7466759693650d-01, &
  0.5334733683389d-08, 0.3594106067745d+01, 0.7489573444450d-01, &
  0.5392665782110d-08, 0.5630254365606d+01, 0.1055449481598d+01, &
  0.6682075673789d-08, 0.1518480041732d+01, 0.2213766559277d+00, &
  0.5079130495960d-08, 0.2739765115711d+01, 0.2132517061319d+00, &
  0.5077759793261d-08, 0.5290711290094d+01, 0.2133464534247d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 151, 160) /&
  0.4832037368310d-08, 0.1404473217200d+01, 0.7160067364790d-01, &
  0.6463279674802d-08, 0.6038381695210d+01, 0.2209183458640d-01, &
  0.6240592771560d-08, 0.1290170653666d+01, 0.3306188016693d+00, &
  0.4672013521493d-08, 0.3261895939677d+01, 0.7796265773310d-01, &
  0.6500650750348d-08, 0.1154522312095d+01, 0.3884652414254d+00, &
  0.6344161389053d-08, 0.6206111545062d+01, 0.7605151500000d-01, &
  0.4682518370646d-08, 0.5409118796685d+01, 0.1073608853559d+01, &
  0.5329460015591d-08, 0.1202985784864d+01, 0.7287631425543d+00, &
  0.5701588675898d-08, 0.4098715257064d+01, 0.8731175355560d-01, &
  0.6030690867211d-08, 0.4132033218460d+00, 0.9846002785331d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 161, 170) /&
  0.4336256312655d-08, 0.1211415991827d+01, 0.4297791515992d+00, &
  0.4688498808975d-08, 0.3765479072409d+01, 0.2127790306879d+00, &
  0.4675578609335d-08, 0.4265540037226d+01, 0.2138191288687d+00, &
  0.4225578112158d-08, 0.5237566010676d+01, 0.3407705765729d+00, &
  0.5139422230028d-08, 0.1507173079513d+01, 0.7233337363710d-01, &
  0.4619995093571d-08, 0.9023957449848d-01, 0.8603097737811d+00, &
  0.4494776255461d-08, 0.5414930552139d+00, 0.7381754420900d-01, &
  0.4274026276788d-08, 0.4145735303659d+01, 0.7574578717200d-01, &
  0.5018141789353d-08, 0.3344408829055d+01, 0.3180992042600d-02, &
  0.4866163952181d-08, 0.3348534657607d+01, 0.7722995774390d-01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 171, 180) /&
  0.4111986020501d-08, 0.4198823597220d+00, 0.1451108196653d+00, &
  0.3356142784950d-08, 0.5609144747180d+01, 0.1274714967946d+00, &
  0.4070575554551d-08, 0.7028411059224d+00, 0.3503323232942d+00, &
  0.3257451857278d-08, 0.5624697983086d+01, 0.5296435984654d+00, &
  0.3256973703026d-08, 0.1857842076707d+01, 0.5297383457582d+00, &
  0.3830771508640d-08, 0.4562887279931d+01, 0.9098186128426d+00, &
  0.3725024005962d-08, 0.2358058692652d+00, 0.1084620721060d+00, &
  0.3136763921756d-08, 0.2049731526845d+01, 0.2346394437820d+00, &
  0.3795147256194d-08, 0.2432356296933d+00, 0.1862120789403d+00, &
  0.2877342229911d-08, 0.5631101279387d+01, 0.1905464808669d+01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 181, 190) /&
  0.3076931798805d-08, 0.1117615737392d+01, 0.3628624111593d+00, &
  0.2734765945273d-08, 0.5899826516955d+01, 0.2131850110243d+00, &
  0.2733405296885d-08, 0.2130562964070d+01, 0.2134131485323d+00, &
  0.2898552353410d-08, 0.3462387048225d+00, 0.5291709230214d+00, &
  0.2893736103681d-08, 0.8534352781543d+00, 0.5302110212022d+00, &
  0.3095717734137d-08, 0.2875061429041d+01, 0.2976424921901d+00, &
  0.2636190425832d-08, 0.2242512846659d+01, 0.1485980103780d+01, &
  0.3645512095537d-08, 0.1354016903958d+01, 0.6044726378023d+00, &
  0.2808173547723d-08, 0.6705114365631d-01, 0.6225157782540d-01, &
  0.2625012866888d-08, 0.4775705748482d+01, 0.5268983110410d-01 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 191, 200) /&
  0.2572233995651d-08, 0.2638924216139d+01, 0.1258454114666d+01, &
  0.2604238824792d-08, 0.4826358927373d+01, 0.2103781122809d+00, &
  0.2596886385239d-08, 0.3200388483118d+01, 0.2162200472757d+00, &
  0.3228057304264d-08, 0.5384848409563d+01, 0.2007689919132d+00, &
  0.2481601798252d-08, 0.5173373487744d+01, 0.1062562936266d+01, &
  0.2745977498864d-08, 0.6250966149853d+01, 0.5651155736444d+00, &
  0.2669878833811d-08, 0.4906001352499d+01, 0.1400015846597d+00, &
  0.3203986611711d-08, 0.5034333010005d+01, 0.7036329877322d+00, &
  0.3354961227212d-08, 0.6108262423137d+01, 0.4549093064213d+00, &
  0.2400407324558d-08, 0.2135399294955d+01, 0.2125476091956d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 201, 210) /&
  0.2379905859802d-08, 0.5893721933961d+01, 0.2140505503610d+00, &
  0.2550844302187d-08, 0.3331940762063d+01, 0.1534957940063d+00, &
  0.2268824211001d-08, 0.1843418461035d+01, 0.2235935264888d+00, &
  0.2464700891204d-08, 0.3029548547230d+01, 0.2091065926078d+00, &
  0.2436814726024d-08, 0.4994717970364d+01, 0.2174915669488d+00, &
  0.2443623894745d-08, 0.2645102591375d+01, 0.1739420156204d+00, &
  0.2318701783838d-08, 0.5700547397897d+01, 0.7530171478090d-01, &
  0.2284448700256d-08, 0.5268898905872d+01, 0.7426161660010d-01, &
  0.2468848123510d-08, 0.5276280575078d+01, 0.2526561439362d+00, &
  0.2814052350303d-08, 0.6130168623475d+01, 0.5636314030725d+00 /
  DATA ((S0(I, J, 2), I = 1, 3), J = 211, NS0Y) /&
  0.2243662755220d-08, 0.6631692457995d+00, 0.8886590321940d-01, &
  0.2330795855941d-08, 0.2499435487702d+01, 0.1056200952181d+01, &
  0.9757679038404d-09, 0.5796846023126d+01, 0.7826370942180d+02 /

  !  SSB-to-Sun, T^1, Y
  DATA ((S1(I, J, 2), I = 1, 3), J = 1, 10) /&
  0.8989047573576d-08, 0.5840593672122d+01, 0.4265981595566d+00, &
  0.7815938401048d-08, 0.1129664707133d+01, 0.2061856251104d+00, &
  0.7550926713280d-08, 0.6196589104845d+00, 0.2204125344462d+00, &
  0.6056556925895d-08, 0.1677494667846d+01, 0.1059381944224d+01, &
  0.5734142698204d-08, 0.4000920852962d+01, 0.5225775174439d+00, &
  0.5614341822459d-08, 0.3486722577328d+01, 0.5368044267797d+00, &
  0.1028678147656d-08, 0.1877141024787d+01, 0.7113454667900d-02, &
  0.7270792075266d-09, 0.5077167301739d+01, 0.6398972393349d+00, &
  0.8734141726040d-09, 0.9069550282609d-01, 0.4194847048887d+00, &
  0.5377371402113d-09, 0.6039381844671d+01, 0.4337116142245d+00 /
  DATA ((S1(I, J, 2), I = 1, 3), J = 11, 20) /&
  0.4729719431571d-09, 0.2153086311760d+01, 0.2132990797783d+00, &
  0.4458052820973d-09, 0.5059830025565d+01, 0.5296909721118d+00, &
  0.4406855467908d-09, 0.2027971692630d+01, 0.1589072916335d+01, &
  0.3101659310977d-09, 0.3317677981860d+01, 0.1052268489556d+01, &
  0.3016749232545d-09, 0.3913703482532d+01, 0.1066495398892d+01, &
  0.3198541352656d-09, 0.1275513098525d+01, 0.1495633313810d+00, &
  0.2142065389871d-09, 0.5301351614597d+01, 0.3163918923335d+00, &
  0.1902615247592d-09, 0.4894943352736d+00, 0.2275259891141d+00, &
  0.1613410990871d-09, 0.2449891130437d+01, 0.1102062672231d+00, &
  0.1576992165097d-09, 0.4211421447633d+01, 0.7626583626240d-01 /
  DATA ((S1(I, J, 2), I = 1, 3), J = 21, 30) /&
  0.1241637259894d-09, 0.4140803368133d+01, 0.5154640627760d+00, &
  0.1313974830355d-09, 0.3591920305503d+01, 0.3664874755930d-01, &
  0.1181697118258d-09, 0.1506314382788d+01, 0.6327837846670d+00, &
  0.1238239742779d-09, 0.7461405378404d+00, 0.3961708870310d-01, &
  0.1010107068241d-09, 0.6271010795475d+00, 0.7329749511860d-01, &
  0.9226316616509d-10, 0.1259158839583d+01, 0.1990721704425d+00, &
  0.8664946419555d-10, 0.3353244696934d+01, 0.5439178814476d+00, &
  0.7757230468978d-10, 0.1447677295196d+01, 0.9491756770005d+00, &
  0.7693168628139d-10, 0.1120509896721d+01, 0.1030928125552d+00, &
  0.5487897454612d-10, 0.4439380426795d+01, 0.8531963191132d+00 /
  DATA ((S1(I, J, 2), I = 1, 3), J = 31, 40) /&
  0.5196118677218d-10, 0.3788856619137d+00, 0.2093666171530d+00, &
  0.5110853339935d-10, 0.1386879372016d+01, 0.2172315424036d+00, &
  0.5027804534813d-10, 0.1647881805466d+00, 0.2164800718209d+00, &
  0.4922485922674d-10, 0.1594315079862d+01, 0.2101180877357d+00, &
  0.6155599524400d-10, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.4447147832161d-10, 0.5480720918976d+01, 0.3235053470014d+00, &
  0.4144691276422d-10, 0.1931371033660d+01, 0.6373574839730d-01, &
  0.4099950625452d-10, 0.5229611294335d+01, 0.6470106940028d+00, &
  0.5060541682953d-10, 0.1731112486298d+01, 0.1422690933580d-01, &
  0.4293615946300d-10, 0.2714571038925d+01, 0.7358765972222d+00 /
  DATA ((S1(I, J, 2), I = 1, 3), J = 41, NS1Y) /&
  0.3545659845763d-10, 0.4451041444634d+01, 0.5265099800692d+00, &
  0.3479112041196d-10, 0.3029385448081d+01, 0.5328719641544d+00, &
  0.3438516493570d-10, 0.2778507143731d+01, 0.8582758298370d-01, &
  0.3297341285033d-10, 0.7898709807584d+00, 0.1104591729320d-01, &
  0.2972585818015d-10, 0.3218785316973d+01, 0.5257585094865d+00, &
  0.2931707295017d-10, 0.4260731012098d+01, 0.5336234347371d+00, &
  0.2897198149403d-10, 0.1120753978101d+01, 0.1173197218910d+00, &
  0.2832293240878d-10, 0.4597682717827d+00, 0.2022531624851d+00, &
  0.2864348326612d-10, 0.2169939928448d+01, 0.9597935788730d-01, &
  0.2852714675471d-10, 0.2377659870578d+01, 0.2118763888447d+01 /

  !  SSB-to-Sun, T^2, Y
  DATA ((S2(I, J, 2), I = 1, 3), J = 1, NS2Y) /&
  0.1609114495091d-11, 0.2831096993481d+01, 0.2061856251104d+00, &
  0.1560330784946d-11, 0.5193058213906d+01, 0.2204125344462d+00, &
  0.1183535479202d-11, 0.5707003443890d+01, 0.5225775174439d+00, &
  0.1158183066182d-11, 0.1782400404928d+01, 0.5368044267797d+00, &
  0.1032868027407d-11, 0.4036925452011d+01, 0.2132990797783d+00, &
  0.6540142847741d-12, 0.4058241056717d+01, 0.4265981595566d+00, &
  0.7305236491596d-12, 0.6175401942957d+00, 0.5296909721118d+00, &
  -0.5580725052968d-12, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.3946122651015d-12, 0.4108265279171d+00, 0.1059381944224d+01 /

  !  SSB-to-Sun, T^0, Z
  DATA ((S0(I, J, 3), I = 1, 3), J = 1, 10) /&
  0.1181255122986d-03, 0.4607918989164d+00, 0.2132990797783d+00, &
  0.1127777651095d-03, 0.4169146331296d+00, 0.5296909721118d+00, &
  0.4777754401806d-04, 0.4582657007130d+01, 0.3813291813120d-01, &
  0.1129354285772d-04, 0.5758735142480d+01, 0.7478166569050d-01, &
  -0.1149543637123d-04, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.3298730512306d-05, 0.5978801994625d+01, 0.4265981595566d+00, &
  0.2733376706079d-05, 0.7665413691040d+00, 0.1059381944224d+01, &
  0.9426389657270d-06, 0.3710201265838d+01, 0.2061856251104d+00, &
  0.8187517749552d-06, 0.3390675605802d+00, 0.2204125344462d+00, &
  0.4080447871819d-06, 0.4552296640088d+00, 0.5225775174439d+00 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 11, 20) /&
  0.3169973017028d-06, 0.3445455899321d+01, 0.5368044267797d+00, &
  0.2438098615549d-06, 0.5664675150648d+01, 0.3664874755930d-01, &
  0.2601897517235d-06, 0.1931894095697d+01, 0.1495633313810d+00, &
  0.2314558080079d-06, 0.3666319115574d+00, 0.3961708870310d-01, &
  0.1962549548002d-06, 0.3167411699020d+01, 0.7626583626240d-01, &
  0.2180518287925d-06, 0.1544420746580d+01, 0.7113454667900d-02, &
  0.1451382442868d-06, 0.1583756740070d+01, 0.1102062672231d+00, &
  0.1358439007389d-06, 0.5239941758280d+01, 0.6398972393349d+00, &
  0.1050585898028d-06, 0.2266958352859d+01, 0.3163918923335d+00, &
  0.1050029870186d-06, 0.2711495250354d+01, 0.4194847048887d+00 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 21, 30) /&
  0.9934920679800d-07, 0.1116208151396d+01, 0.1589072916335d+01, &
  0.1048395331560d-06, 0.3408619600206d+01, 0.1021328554739d+02, &
  0.8370147196668d-07, 0.3810459401087d+01, 0.2535050500000d-01, &
  0.7989856510998d-07, 0.3769910473647d+01, 0.7329749511860d-01, &
  0.5441221655233d-07, 0.2416994903374d+01, 0.1030928125552d+00, &
  0.4610812906784d-07, 0.5858503336994d+01, 0.4337116142245d+00, &
  0.3923022803444d-07, 0.3354170010125d+00, 0.1484170571900d-02, &
  0.2610725582128d-07, 0.5410600646324d+01, 0.6327837846670d+00, &
  0.2455279767721d-07, 0.6120216681403d+01, 0.1162474756779d+01, &
  0.2375530706525d-07, 0.6055443426143d+01, 0.1052268489556d+01 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 31, 40) /&
  0.1782967577553d-07, 0.3146108708004d+01, 0.8460828644453d+00, &
  0.1581687095238d-07, 0.6255496089819d+00, 0.3340612434717d+01, &
  0.1594657672461d-07, 0.3782604300261d+01, 0.1066495398892d+01, &
  0.1563448615040d-07, 0.1997775733196d+01, 0.2022531624851d+00, &
  0.1463624258525d-07, 0.1736316792088d+00, 0.3516457698740d-01, &
  0.1331585056673d-07, 0.4331941830747d+01, 0.9491756770005d+00, &
  0.1130634557637d-07, 0.6152017751825d+01, 0.2968341143800d-02, &
  0.1028949607145d-07, 0.2101792614637d+00, 0.2275259891141d+00, &
  0.1024074971618d-07, 0.4071833211074d+01, 0.5070101000000d-01, &
  0.8826956060303d-08, 0.4861633688145d+00, 0.2093666171530d+00 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 41, 50) /&
  0.8572230171541d-08, 0.5268190724302d+01, 0.4110125927500d-01, &
  0.7649332643544d-08, 0.5134543417106d+01, 0.2608790314060d+02, &
  0.8581673291033d-08, 0.2920218146681d+01, 0.1480791608091d+00, &
  0.8430589300938d-08, 0.3604576619108d+01, 0.2172315424036d+00, &
  0.7776165501012d-08, 0.3772942249792d+01, 0.6373574839730d-01, &
  0.8311070234408d-08, 0.6200412329888d+01, 0.3235053470014d+00, &
  0.6927365212582d-08, 0.4543353113437d+01, 0.8531963191132d+00, &
  0.6791574208598d-08, 0.2882188406238d+01, 0.7181332454670d-01, &
  0.5593100811839d-08, 0.1776646892780d+01, 0.7429900518901d+00, &
  0.4553381853021d-08, 0.3949617611240d+01, 0.7775000683430d-01 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 51, 60) /&
  0.5758000450068d-08, 0.3859251775075d+01, 0.1990721704425d+00, &
  0.4281283457133d-08, 0.1466294631206d+01, 0.2118763888447d+01, &
  0.4206935661097d-08, 0.5421776011706d+01, 0.1104591729320d-01, &
  0.4213751641837d-08, 0.3412048993322d+01, 0.2243449970715d+00, &
  0.5310506239878d-08, 0.5421641370995d+00, 0.5154640627760d+00, &
  0.3827450341320d-08, 0.8887314524995d+00, 0.1510475019529d+00, &
  0.4292435241187d-08, 0.1405043757194d+01, 0.1422690933580d-01, &
  0.3189780702289d-08, 0.1060049293445d+01, 0.1173197218910d+00, &
  0.3226611928069d-08, 0.6270858897442d+01, 0.2164800718209d+00, &
  0.2893897608830d-08, 0.5117563223301d+01, 0.6470106940028d+00 /
  DATA ((S0(I, J, 3), I = 1, 3), J = 61, NS0Z) /&
  0.3239852024578d-08, 0.4079092237983d+01, 0.2101180877357d+00, &
  0.2956892222200d-08, 0.1594917021704d+01, 0.3092784376656d+00, &
  0.2980177912437d-08, 0.5258787667564d+01, 0.4155522422634d+00, &
  0.3163725690776d-08, 0.3854589225479d+01, 0.8582758298370d-01, &
  0.2662262399118d-08, 0.3561326430187d+01, 0.5257585094865d+00, &
  0.2766689135729d-08, 0.3180732086830d+00, 0.1385174140878d+00, &
  0.2411600278464d-08, 0.3324798335058d+01, 0.5439178814476d+00, &
  0.2483527695131d-08, 0.4169069291947d+00, 0.5336234347371d+00, &
  0.7788777276590d-09, 0.1900569908215d+01, 0.5217580628120d+02 /

  !  SSB-to-Sun, T^1, Z
  DATA ((S1(I, J, 3), I = 1, 3), J = 1, 10) /&
  0.5444220475678d-08, 0.1803825509310d+01, 0.2132990797783d+00, &
  0.3883412695596d-08, 0.4668616389392d+01, 0.5296909721118d+00, &
  0.1334341434551d-08, 0.0000000000000d+00, 0.0000000000000d+00, &
  0.3730001266883d-09, 0.5401405918943d+01, 0.2061856251104d+00, &
  0.2894929197956d-09, 0.4932415609852d+01, 0.2204125344462d+00, &
  0.2857950357701d-09, 0.3154625362131d+01, 0.7478166569050d-01, &
  0.2499226432292d-09, 0.3657486128988d+01, 0.4265981595566d+00, &
  0.1937705443593d-09, 0.5740434679002d+01, 0.1059381944224d+01, &
  0.1374894396320d-09, 0.1712857366891d+01, 0.5368044267797d+00, &
  0.1217248678408d-09, 0.2312090870932d+01, 0.5225775174439d+00 /
  DATA ((S1(I, J, 3), I = 1, 3), J = 11, NS1Z) /&
  0.7961052740870d-10, 0.5283368554163d+01, 0.3813291813120d-01, &
  0.4979225949689d-10, 0.4298290471860d+01, 0.4194847048887d+00, &
  0.4388552286597d-10, 0.6145515047406d+01, 0.7113454667900d-02, &
  0.2586835212560d-10, 0.3019448001809d+01, 0.6398972393349d+00 /

  !  SSB-to-Sun, T^2, Z
  DATA ((S2(I, J, 3), I = 1, 3), J = 1, NS2Z) /&
  0.3749920358054d-12, 0.3230285558668d+01, 0.2132990797783d+00, &
  0.2735037220939d-12, 0.6154322683046d+01, 0.5296909721118d+00 /

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   t = (xjd-xjd2000)/yearj
   t2= t * t
   do j = 1, 3 ! loop on x, y, z
      aa = 0d0
      bb = 0d0

      do k = 1, ns0(j)
         arg = s0(2, k, j) + s0(3, k, j) * t
         amp = s0(1, k, j)
         aa = aa + amp * cos(arg)
         bb = bb - s0(3, k, j) * amp * sin(arg)
      enddo

      do k = 1, ns1(j)
         arg = s1(2, k, j) + s1(3, k, j) * t
         cc = cos(arg)
         amp = s1(1, k, j)
         aa = aa + amp * t * cc
         bb = bb + (amp * cc - s1(3, k, j) * amp * t * sin(arg))
      enddo

      do k = 1, ns2(j)
         arg = s2(2, k, j) + s2(3, k, j) * t
         cc = cos(arg)
         amp = s2(1, k, j)
         aa = aa + amp * t2 * cc
         bb = bb + (2d0 * t * amp * cc - s2(3, k, j) * amp * t2 * sin(arg))
      enddo

      xyz(j) = aa ! position in au
      xyzp(j) = bb/yearj ! velocity in au/day
   enddo
   return
end subroutine baryc_h


  !****************************************************************************** 
  subroutine planet_helio_simon(iref, iplan, datjd, xlong, xlat, dist, pos, vit)
  !******************************************************************************     
  !
  !      Heliocentric cooordinates of the planets Mercury to Neptune and EMB
  !      and geocentric for the Sun
  !      Ecliptic J2000 (iref = 0) or mean ecliptic and equinox of date (iref =1)
  !      medium precision over the interval  500- 2500
  !      Heliocentric geometric coordinates (geocentric for the Sun)
  !
  !      angular units  : degrees
  !
  !      Author  : F. Mignard  September 2015 from core routine provided by J.L. Simon
  !
  !      Source  : 
  !      The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront, &
  !      M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
  !      Longitudes, Paris, France).  From comparisons with JPL
  !      ephemeris DE102, they quote the following maximum errors
  !      over the interval 1800-2050:
  !
  !                     L (arcsec)    B (arcsec)      R (km)
  !
  !        Mercury          4             1             300
  !        Venus            5             1             800
  !        EMB              6             1            1000
  !        Mars            17             1            7700
  !        Jupiter         71             5           76000
  !        Saturn          81            13          267000
  !        Uranus          86             7          712000
  !        Neptune         11             1          253000
  !
  !     Over the interval 1000-3000, they report that the accuracy is no
  !     worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
  !     accuracy declines.
  !
  !******************************************************
  !
  !***  INPUT
  !
  !      iref        : = 0 Ecliptic J2000 ;  =1 : mean ecliptic of date and rotating frame for the velocity
  !      iplan       : planet id (0 : Sun (from Earth), 1 : Mercury ... 8: Neptune, 9 :EMB)
  !      datjd       : date in julian days
  !
  !***  OUTPUT
  !
  !      xlong       : ecliptic longitude in degrees
  !      xlat        : ecliptic latitude  in degrees
  !      dist        : distance to the sun in au
  !      pos         : position vector pos(3) in au
  !      vit         : velocity vector vit(3) in au/day
  !
  !******************************************************


  implicit none

  integer, parameter :: iframe = 1 ! ecliptic frame for precession and barycentric orbit
  integer, parameter :: iearth = 3 !


  integer, intent(in)           :: iref, iplan
  real(kind = dp), intent(in)   :: datjd
  real(kind = dp), intent(out)  :: xlong, xlat, dist, pos(3), vit(3)

  real(kind = dp), dimension(3) :: dvit, rr, vv, xyz, xyzp

  
  if(iplan==0) then  ! sun  from the geocenter
     call ephem_planet_simon(datjd, iearth, iframe, rr, vv)  ! ICRS or ecplitic J2000, km, km/s
     call baryc_m(datjd, xyz, xyzp)                         ! ecliptic J2000 km, km/s
     pos = -(rr/xau - xyz)
     vit = -(vv/audtokms - xyzp)
  else   
     call ephem_planet_simon(datjd, iplan, iframe, rr, vv)  ! ICRS or ecplitic J2000, km, km/s
     call baryc_m(datjd, xyz, xyzp)                         ! ecliptic J2000 au au/days
     pos = rr/xau - xyz
     vit = vv/audtokms - xyzp
  endif  
  !
  ! Position and velocity vectors at xjd : mean ecliptic of date
  ! Warning : velocity transformation accounts for the rotating frame dragging velocity
  !

  if (iref == 1) then ! mean ecliptic of date
     call preces(iframe, xjd2000, datjd, pos, pos) !pos in au
     call preces(iframe, xjd2000, datjd, vit, vit) !vit in au/day
     call v_preces(iframe, xjd2000, datjd, pos, dvit) !dvit in au/day :: rotational velocity from rotating frame
     vit = vit + dvit
  endif
  call carsphe(pos, dist, xlong, xlat)
  

  return
  end subroutine planet_helio_simon   
  
!!  Utility routines takedn from FM astro and math modules and used in this module. This allows it to be self contained  
!
!************************************************************
 function sind(x)
 !***********************************************************
      !
      !     sin with argument in degree
      !     could be deleted if sind is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, sind
      sind = sin(x * degrad)
   end function sind

   !***********************************************************
   function cosd(x)
      !***********************************************************
      !
      !     cos with argument in degree
      !     could be deleted if cosd is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, cosd
      cosd = cos(x * degrad)
   end function cosd


   !***********************************************************
   function tand(x)
      !***********************************************************
      !
      !     tan with argument in degree
      !     could be deleted if tand is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, tand
      tand = tan(x * degrad)
   end function tand

   !***********************************************************
   function asind(x)
      !***********************************************************
      !
      !     asin with result in degree
      !     could be deleted if asind is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, asind
      asind = asin(x) * raddeg
   end function asind
   !***********************************************************
   function acosd(x)
      !***********************************************************
      !
      !     acos with result in degree
      !     could be deleted if acosd is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, acosd
      acosd = acos(x) * raddeg
   end function acosd

   !***********************************************************
   function atand(x)
      !***********************************************************
      !
      !     atan with result in degree
      !     could be deleted if atand is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, atand
      atand = atan(x) * raddeg
   end function atand


   !***********************************************************
   function atan2d(x, y)
      !***********************************************************
      !
      !     atan2 with result in degree
      !     could be deleted if atan2d is accepted by the Fortran implementation
      !
      implicit none
      real(kind = dp) x, y, atan2d
      atan2d = atan2(x, y) * raddeg
   end function atan2d
   !******************************************************
   function deg(xrad)
      !******************************************************
      !
      !
      !       Transformation from radians to degrees
      !
      !       author F. Mignard OCA/CERGA
      !
      !*** Input
      !
      !       xrad  : angle in radians
      !
      !*** Output
      !
      !       deg(xrad) : angle in degrees
      !
      !
      !******************************************************
      implicit none

      real(kind = dp), intent(in) :: xrad
      real(kind = dp) :: deg

      real(kind = dp), parameter :: pi = 3.141592653589793238462643d0

      deg = (180.0d0/pi) * xrad

      return
   end function deg
   !************************************************************************!******************************************************
   function rad(xdeg)
      !******************************************************
      !
      !
      !       Transformation from degree to radians
      !
      !       author F. Mignard OCA/CERGA
      !
      !*** Input
      !
      !       xdeg  : angle in decimal degrees
      !
      !*** Output
      !
      !       rad(xdeg) : angle in radians
      !
      !
      !******************************************************
      implicit none
      real(kind = dp), intent(in) :: xdeg
      real(kind = dp) :: rad

      real(dp), parameter :: pi = 3.141592653589793238462643d0

      rad = pi/180.0d0 * xdeg

      return
   end function rad   
   !******************************************************
   function sqr(x)
      !******************************************************
      !       square of a real number (for compatibility with Pascal source)
      !
      !       author F. Mignard OCA/CERGA
      !
      !*** Input
      !
      !       x  : number to be squared
      !
      !*** Output
      !
      !      sqr(x) : x*x
      !
      !
      !******************************************************
      implicit none
      real(kind = dp) :: x, sqr

      sqr = x * x

      return
   end function sqr  
   !******************************************************
   function angle(x, y)
   !******************************************************
      !
      !
      !     angle solution of x = r*cos(teta), y = r*sin(teta)  r >0
      !
      !     author F. Mignard OCA/CERGA
      !
      !***  INPUT
      !
      !     x      : r*cos(teta)
      !     y      : r*sin(teta)
      !
      !***  OUTPUT
      !
      !     angle  : angle in degrees in 0 - 360
      !
      !******************************************************
      implicit none

      real(kind = dp), intent(in) :: x, y
      real(kind = dp) :: angle
      real(kind = dp) :: xx

      if ((x * x + y * y) == 0) then
         xx = 0 ! conventional value for the zero vector
      else
         xx = deg(atan2(y, x)) ! atan2d is not in double precision
      endif

      if (xx < 0) then
         xx = xx + 360.d0
      endif

      angle = xx
      return
   end function angle   
  !******************************************************
   function fp(x)
      !******************************************************
      !
      !     fp    : function fractional part of x
      !     0 <= fp(x)  <1   with x > 0 or  < 0
      !
      !     f. mignard  december 1998
      !
      !     fp( 3.14)  = 0.14
      !     fp(-3.14)  = 0.86
      !     fp( 2.71)  = 0.71
      !     fp(-2.71)  = 0.29
      !******************************************************
      implicit none

      real(kind = dp) :: fp
      real(kind = dp), intent(in) :: x
      !
      if ((x - int(x)) .eq. 0) then
         fp = 0.
      else
         fp = x - floor(x)
      endif
      return
   end function fp  
   !*****************************************************************
   function xmodpi(x, y)
      !*****************************************************************
      !
      !     Modulo of x within [-y/2, +y/2]
      !
      !*****************************************************************
      !    INPUT
      !          x            number wose remainder is searched
      !          y            module
      !
      !    OUTPUT
      !          xmodpi       remainder in   [-y/2, +y/2]
      !
      !*****************************************************************
      implicit none
      real(kind = dp) :: xmodpi
      real(kind = dp), intent(in) :: x, y

      real(kind = dp) :: r


      r = modulo(x, y)

      if (r > 0.5 * y) then
         r = r - y
      endif

      xmodpi = r
      return
   end function xmodpi  
!******************************************************
subroutine carsphe(v, rho, alpha, delta)
   !******************************************************
   !
   !     Transformation from cartesian coordinates to spherical.
   !     latitude from the equator and not from the pole
   !
   !     author F. Mignard OCA/CERGA
   !
   !***  INPUT
   !
   !       v    : input vector v(3)
   !
   !***  OUTPUT
   !
   !      rho   : modulus of v
   !      alpha : longitude in degrees
   !      delta : latitude in degrees
   !
   !      for the poles : alpha = 0
   !******************************************************
   implicit none

   real(kind = dp) v(3), rho, alpha, delta, rhop

   rho = sqrt(dot_product(v, v))
   rhop = sqrt(v(1) * v(1) + v(2) * v(2))

   if (rho == 0) then
      alpha = 0d0 ! conventional values
      delta = 0d0
   else
      if (rhop == 0) then
      if (v(3) > 0) then
         delta = 90.0d0
      else
         delta = -90.0d0
      endif
      alpha = 0.0d0 ! conventional value
   else
      delta = deg(atan(v(3)/rhop))
      alpha = angle(v(1), v(2))
      endif
   endif
   return
end subroutine carsphe
!******************************************************
subroutine sphecar(rho, alpha, delta, v)
   !******************************************************
   !
   !     Transformation from  spherical to cartesian coordinates
   !     latitude from the equator and not from the pole
   !
   !     author F. Mignard OCA/CERGA
   !
   !***  INPUT
   !
   !      rho   : radius vector
   !      alpha : longitude in degrees
   !      delta : latitude in degrees
   !
   !***  OUTPUT
   !
   !       v    : input vector v(3)
   !
   !******************************************************
   implicit none

   real(kind = dp) v(3), rho, alpha, delta, cd

   cd = cosd(delta)
   v(1) = rho * cosd(alpha) * cd
   v(2) = rho * sind(alpha) * cd
   v(3) = rho * sind(delta)

   return
end subroutine sphecar
!******************************************************
function unitvect(v)
   !******************************************************
   !
   !     UNit vector of a 3-D vector v
   !
   ! author  F. Mignard OCA/Lagrange
   !
   !***  INPUT
   !
   !     v    : input vector v(3)
   !
   !***  OUTPUT
   !
   !     unitvect     : v/sqrt(v.v)
   !******************************************************
   implicit none

   real(kind = dp),intent(in)   :: v(3)
   real(kind = dp)              :: unitvect(3)
   real(kind = dp)              :: xnorm

   xnorm = dot_product(v,v)

   if(xnorm > 0d0) then
      unitvect = v/sqrt(dot_product(v,v))
   else
      print * , 'vector with zero norm found in function unitvect'
      print *,  'program is aborted '
      stop
   endif

   return
end function unitvect
!******************************************************
function normvect(v)
   !******************************************************
   !
   !    Euclidan norm of  a 3-D vector
   !
   ! author  F. Mignard OCA/Lagrange
   !
   !***  INPUT
   !
   !     v    : input vector v(3)
   !
   !***  OUTPUT
   !
   !     normvect     :sqrt(v.v)
   !******************************************************
   implicit none

   real(kind = dp),intent(in)   :: v(:)
   real(kind = dp)              :: normvect

   normvect = sqrt(dot_product(v,v))

   return
end function normvect

!******************************************************
subroutine matrot(angle, iaxe, xmrot)
   !***************************************************
   !
   !
   !    Canonical rotation matrix around one axis
   !
   !
   ! Author      : F. Mignard OCA/Cassiopee
   !
   !*** INPUT
   !
   !    angle    : angle of rotation in degrees
   !    iaxe     : 1 (=x) , 2 (=y) , 3(=z) for the axis of rotation
   !
   !***OUTPUT
   !
   !    xmrot    : rotation matrix 3x3
   !
   !
   !  convention : the reference frame is rotated, not the vector:
   !               so by applying the matrix to a vector in the initial frame
   !               one gets the coordinate of the same vector in the rotated frame.
   !
   !******************************************************
   implicit none

   real(kind = dp), intent(in) :: angle
   integer, intent(in) :: iaxe

   real(kind = dp), dimension(3, 3), intent(out) :: xmrot

   real(kind = dp) :: cc, ss

   cc = cosd(angle)
   ss = sind(angle)
   xmrot = 0.0d0

   if ((iaxe > 3).or.(iaxe < 1)) then
      write(*, *) 'axis value in rota : ', iaxe
      write(*, *) 'execution terminated'
      stop
   else
      select case (iaxe)
      case(1)
         xmrot(1, 1) = 1d0
         xmrot(2, 2) = cc; xmrot(2, 3) = ss
         xmrot(3, 2) = -ss; xmrot(3, 3) = cc
      case(2)
         xmrot(1, 1) = cc; xmrot(1, 3) = -ss
         xmrot(2, 2) = 1d0
         xmrot(3, 1) = ss; xmrot(3, 3) = cc
      case(3)
         xmrot(1, 1) = cc; xmrot(1, 2) = ss
         xmrot(2, 1) = -ss; xmrot(2, 2) = cc
         xmrot(3, 3) = 1d0
      endselect
   endif
   return
end subroutine matrot
   !******************************************************
   subroutine rota(v_in, angle, iaxe, v_out)
      !******************************************************
      !
      !
      !    Application of a canonical rotation to a vector
      !
      !
      ! Author      : F. Mignard OCA/CERGA
      !
      !*** INPUT
      !
      !    v_in     : input 3D vector  array(3)
      !    angle    : angle of rotation in degrees
      !    iaxe     : 1 (=x) , 2 (=y) , 3(=z) for the axis of rotation
      !
      !***OUTPUT
      !
      !    v_out     : output 3D vector  array(3)
      !
      !
      !    remark   : the routine may be called with the same input and output vector like
      !              call rota(v,angle, iaxe, v).
      !
      !  convention : the reference frame is rotated, not the vector  : so one has
      !               the coordinate of the same vector in the old and new frame.
      !
      !******************************************************
      implicit none

      real(kind = dp), dimension(3), intent(in) :: v_in
      real(kind = dp), intent(in) :: angle
      integer, intent(in) :: iaxe

      real(kind = dp), dimension(3), intent(out) :: v_out

      real(kind = dp), dimension(3) :: v
      real(kind = dp) :: cc, ss

      !********************************************************

      cc = cosd(angle)
      ss = sind(angle)
      v = v_in ! allow to call :  call rota(V,angle, iaxe, V)
      if ((iaxe > 3).or.(iaxe < 1)) then
         write(*, *) 'axis value in rota : ', iaxe
         write(*, *) 'execution terminated'
         stop
      else
         select case (iaxe)
         case(1)
            v_out(1) = v(1)
            v_out(2) = cc * v(2) + ss * v(3)
            v_out(3) = -ss * v(2) + cc * v(3)
         case(2)
            v_out(1) = cc * v(1) - ss * v(3)
            v_out(2) = v(2)
            v_out(3) = ss * v(1) + cc * v(3)
         case(3)
            v_out(1) = cc * v(1) + ss * v(2)
            v_out(2) = -ss * v(1) + cc * v(2)
            v_out(3) = v(3)
         endselect
      endif
      return
   end subroutine rota
!******************************************************
subroutine left_rot(xmatin, angle, iaxe, xmatout)
   !***************************************************
   !
   !
   !    Rotation of a matrix with a left canonical rotation matrix
   !
   !
   ! Author      : F. Mignard OCA/Cassiopee
   !
   !*** INPUT
   !
   !    xmatin   : 3x3 rotation matrix in input
   !    angle    : angle of rotation in degrees
   !    iaxe     : 1 (=x) , 2 (=y) , 3(=z) for the axis of rotation
   !
   !***OUTPUT
   !
   !    xmatout  : rotated matrix 3x3   xmatout : R(iaxe,angle)*xmatin
   !
   !
   !  convention : the reference frame is rotated, not the vector:
   !               so by applying the matrix to a vector in the initial frame
   !               one gets the coordinate of the same vector in the rotated frame.
   !
   !******************************************************
implicit none

real(kind = dp), dimension(3, 3), intent(in) :: xmatin
real(kind = dp), intent(in) :: angle
integer, intent(in) :: iaxe

real(kind = dp), dimension(3, 3), intent(out):: xmatout

real(kind = dp), dimension(3, 3) :: ymat , zmat


zmat = xmatin    ! allow :  call left_rot(xmat,angle, iaxe, xmat)

if ((iaxe > 3).or.(iaxe < 1)) then
   write(*, *) 'axis value in left_rot : ', iaxe
   write(*, *) 'execution terminated'
   stop
else
   call matrot(angle, iaxe, ymat)
   xmatout = matmul(ymat, zmat)
endif
return
end subroutine left_rot

   !*****************************************************************
   subroutine precession_equ_2006(t, xmatrot)
      !*****************************************************************
      !    Percession iau 2006
      !    version 1 :  September 2015
      !    source : Table 4.2.1 in CDT Explanation
      !    computation of the precession matrix for transformation between
      !    the initial epoch (J2000) and a final (t) in equatorial coordinates
      !
      !
      !    F. Mignard for IAU 2006 Precession model
      !
      !    INPUT
      !
      !      t       :  date in julian days
      !
      !    OUTPUT
      !
      !     xmatrot  : matrix(3,3) between the two reference frames.
      !
      !                convention   for the application :
      !                v(t) = xmatrot*v(J2000)
      !*****************************************************************
      !

      implicit none

      real(kind = dp), intent(in) :: t
      real(kind = dp), dimension(3, 3), intent(out) :: xmatrot

      real(kind=dp), parameter, dimension(3,3) :: xmident = reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0] , shape(xmident))
      real(kind=dp), parameter :: tiny = 1d-10


      real(kind = dp) :: tt, dzeta, za, teta

       tt = (t - xjd2000) / century

      !
      !        values of the three precession angles  IAU 2006
      !

      if(abs(tt)> tiny) then  ! if the call time is J2000 ==> identify matrix is returned without computation
          dzeta =     (+2.650545d0 + tt*(2306.083227d0  +tt*( 29.88499d-2 + tt*( 18.01828d-3  +tt*(-0.05971d-4 +tt*(-0.03173d-5))))))*asdeg
          za    =     (-2.650545d0 + tt*(2306.077181d0  +tt*(109.27348d-2 + tt*( 18.26837d-3  +tt*(-0.28596d-4 +tt*(-0.02904d-5))))))*asdeg
          teta  =     (0d0         + tt*(2004.191903d0  +tt*(-42.94934d-2 + tt*(-41.82264d-3  +tt*(-0.07089d-4 +tt*(-0.01274d-5))))))*asdeg

          call matrot(-dzeta, 3, xmatrot)
          call left_rot(xmatrot, teta, 2, xmatrot)
          call left_rot(xmatrot,  -za, 3, xmatrot)
       else
          xmatrot = xmident  ! identity matrix(3,3)
       endif

      return
   end subroutine precession_equ_2006
  subroutine v_preces(iframe, ti, tf, vi, dv)
      !*****************************************************************
      !
      !Rotating velocity caused by the precession
      !Given a position vector vi in the frame at t1 (e.g equinox of date) the
      !program returns the correction of the velocity vector omega X vi  of the particle
      !due to the rotating frame.
      !dv is given in the frame used at tf (e.g J2000 frame) but vi must be in the initial frame.
      !
      ! This is useful to determine the velocity vector of a planet in the J2000 frame
      ! from an ephemeris given in the mean equinox of date.
      ! V1 with IAU 76 Precession
      ! V2 with IAU 2006 precession September 2015
      !
      ! F. Mignard OCA/Lagrange 2005, 2015
      !
      !iframe = 0  : equatorial iframe  1 : ecliptic frame
      !ti : initial time in julian days  corresponding to vi
      !tf : final time in julian days corresponding to vf
      !
      !    INPUT
      !
      !      iframe  : integer 0 : equator, 1 : ecliptic
      !      ti      : initial date in julian days
      !      tf      : final   date in julian days
      !      vi      : initial vector in the frame at ti
      !
      !    OUTPUT
      !
      !      dv      : velocity due to the rotation of the frame in unit of vi per day.
      !
      !*****************************************************************
      !      use constantes
      implicit none

      real(kind = dp), parameter :: step = 5d0 * yearj ! step for the numerical derivative of precession in days

      integer, intent(in) :: iframe
      real(kind = dp), intent(in) :: ti, tf
      real(kind = dp), dimension(3), intent(in) :: vi
      real(kind = dp), dimension(3), intent(out) :: dv

      real(kind = dp), dimension(3, 3) :: dotmat
      real(kind = dp), dimension(3, 3) :: xmatrot1, xmatrot2, xmatrot3, xmatrot4
      real(kind = dp), dimension(3) :: v
      integer :: i, j


      v = vi ! allow to call :  call preces(iframe,ti,tf, v, v)

      !
      ! Derivative of the precession matrix with respect to tf
      !

      if(iframe == 0) then

      call mat_prec_equ(ti, tf - step, xmatrot1)
      call mat_prec_equ(ti, tf + step, xmatrot2)
      call mat_prec_equ(ti, tf - 2d0 * step, xmatrot3)
      call mat_prec_equ(ti, tf + 2d0 * step, xmatrot4)

      else if (iframe == 1) then

         call mat_prec_ecl(ti, tf - step, xmatrot1)
         call mat_prec_ecl(ti, tf + step, xmatrot2)
         call mat_prec_ecl(ti, tf - 2d0 * step, xmatrot3)
         call mat_prec_ecl(ti, tf + 2d0 * step, xmatrot4)

      else
         print *, 'program interrupted in procedure v_preces  with iframe = ', iframe
         print *, 'should be 0 (equatorial) or  1  (ecliptic)'
         stop
      endif

      dotmat = ((xmatrot3 - xmatrot4) + 8d0 * (xmatrot2 - xmatrot1))/(12d0 * step) ! derivative in radians per day

      dv(1) = 0d0
      dv(2) = 0d0
      dv(3) = 0d0

      do i = 1, 3
         do j = 1, 3
            dv(i) = dv(i) + dotmat(i, j) * v(j)
         enddo
      enddo
      return
   end subroutine v_preces
   
   !*****************************************************************
   subroutine precession_ecl_2006(t, xmatrot)
      !*****************************************************************
      !    Prcession iau 2006
      !    version 1 :  September 2015
      !    source : Table 4.2.1 in CDT Explanation
      !    computation of the precession matrix for transformation between
      !    the initial epoch (J2000) and a final (t) in ecliptic coordinates
      !
      !    F. Mignard for IAU 2006 Precession model
      !
      !    INPUT
      !
      !      t       :  date in julian days
      !
      !    OUTPUT
      !
      !     xmatrot  : matrix(3,3) between the two reference frames.
      !
      !                convention   for the application :
      !                v(t) = xmatrot*v(J2000)
      !*****************************************************************
      !

      implicit none

      real(kind = dp), intent(in) :: t
      real(kind = dp), dimension(3, 3), intent(out) :: xmatrot

      real(kind=dp), parameter, dimension(3,3) :: xmident = reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0] , shape(xmident))
      real(kind=dp), parameter :: tiny = 1d-10


      real(kind = dp) :: tt, gpia, pia, pa
 !     real(kind = dp) :: c1, c2, c3, s1, s2, s3

       tt = (t - xjd2000) / century

      !
      !        values of the three precession angles  IAU 2006
      !

      if(abs(tt)> tiny) then  ! if the call time is J2000 ==> identify matrix is returned without computation
         gpia = 174.874109333d0 + tt*(-867.95758d0 +tt*(0.157992d0  +tt*(-0.5371d-3 +tt*(-0.4797d-4  +tt*0.0072d-5))))*asdeg
         pia  =                   tt*(46.998973d0 +tt*(-3.34926d-2 +tt*(-0.12559d-3 +tt*(0.00113d-4 +tt*(-0.00022d-5)))))*asdeg
         pa  =                    tt*(5028.796195d0 +tt*(110.54348d-2 +tt*(0.07964d-3 +tt*(-0.23857d-4 +tt*(-0.00383d-5)))))*asdeg


         call matrot(gpia, 3, xmatrot)
         call left_rot(xmatrot, pia, 1, xmatrot)
         call left_rot(xmatrot, -gpia-pa, 3, xmatrot)
       else
          xmatrot = xmident   ! identity matrix(3,3)
       endif
!       call mat_ecl_old(xjd2000,t , xmatrot)  ! test with old routine
      return
   end subroutine precession_ecl_2006
     !*****************************************************************
   subroutine mat_prec_ecl(ti, tf, xmat)
      !*****************************************************************
      !    computation of the precession matrix for transformation between
      !    an initial epoch (ti) and a final (tf) in ecliptic coordinates
      !    F. Mignard  September 2015 for this version with IAU 2006
      !
      !    INPUT
      !
      !      ti      : initial date in julian days
      !      tf      : final   date in julian days
      !
      !    OUTPUT
      !
      !      xmat    : matrix(3,3) between the two reference frames.
      !                mean ecliptic at ti and mean ecliptic at tf
      !                convention   for the application :
      !                v(tf) = a*v(ti)
      !*****************************************************************

   implicit none

   real(kind = dp), intent(in) :: ti, tf
   real(kind = dp), dimension(3, 3), intent(out) :: xmat

   real(kind = dp), dimension(3, 3)              :: xmat1, xmat2

   call precession_ecl_2006(ti, xmat1)  ! no computation if ti = J2000 ID matrix returned

   call precession_ecl_2006(tf, xmat2)  ! no computation if tf = J2000 ID matrix returned

   xmat = matmul(xmat2, transpose(xmat1)) ! from ti to J2000, then from J2000 to tf

   return
   end subroutine mat_prec_ecl
     !*****************************************************************
   subroutine mat_prec_equ(ti, tf, xmat)
      !*****************************************************************
      !    computation of the precession matrix for transformation between
      !    an initial epoch (ti) and a final (tf) in equatorial coordinates
      !    F. Mignard  September 2015 for this version with IAU 2006
      !
      !    INPUT
      !
      !      ti      : initial date in julian days
      !      tf      : final   date in julian days
      !
      !    OUTPUT
      !
      !      xmat    : matrix(3,3) between the two reference frames.
      !                mean ecliptic at ti and mean ecliptic at tf
      !                convention   for the application :
      !                v(tf) = a*v(ti)
      !*****************************************************************

   implicit none

   real(kind = dp), intent(in) :: ti, tf
   real(kind = dp), dimension(3, 3), intent(out) :: xmat

   real(kind = dp), dimension(3, 3)              :: xmat1, xmat2

   call precession_equ_2006(ti, xmat1) ! no computation if ti = J2000 ID matrix returned

   call precession_equ_2006(tf, xmat2) ! no computation if tf = J2000 ID matrix returned

   xmat = matmul(xmat2, transpose(xmat1)) ! from ti to J2000, then from J2000 to tf

   return
   end subroutine mat_prec_equ
!*****************************************************************
   subroutine preces(iframe, ti, tf, vi, vf)
!*****************************************************************
     !
      !Change of reference frame from precession between ti and tf
      !iframe = 0  : equatorial iframe  1 : ecliptic frame
      !ti : initial time in julian days  corresponding to vi
      !tf : final time in julian days corresponding to vf
      ! V1 with IAU 76 Precession
      ! V2 with IAU 2006 precession September 2015
	  ! F. Mignard OCA/Cassiopee  2004, 2015
      !
      !    INPUT
      !
      !      iframe  : integer 0 : equator, 1 : ecliptic
      !      ti      : initial date in julian days
      !      tf      : final   date in julian days
      !      vi      : initial vector in the frame at ti
      !
      !    OUTPUT
      !
      !      vf      : final  vector in the frame at tf
      !
      !*****************************************************************
      implicit none

      integer, intent(in) :: iframe
      real(kind = dp), intent(in) :: ti, tf
      real(kind = dp), dimension(3), intent(in) :: vi
      real(kind = dp), dimension(3), intent(out) :: vf

      real(kind = dp), dimension(3, 3) :: xmatrot
      real(kind = dp), dimension(3) :: v
      integer :: i, j


      v = vi ! allow to call :  call preces(iframe,ti,tf, v, v)

      if (iframe == 0) then
         call mat_prec_equ(ti, tf, xmatrot)
      else if (iframe == 1) then
         call mat_prec_ecl(ti, tf, xmatrot)
      else
         print *, 'program interrupted in procedure preces  with iframe = ', iframe
         stop
      endif

      do i = 1, 3
         vf(i) = 0d0
         do j = 1, 3
            vf(i) = vf(i) + xmatrot(i, j) * v(j)
         enddo
      enddo
      return
   end subroutine preces
   !******************************************************
   function obliquity(datjd)
      !******************************************************
      !
      !     Obliquity of the ecliptic at the JD date datjd
      !     unit  : degrees
      !     F. Mignard  August 2001
      !
      !     INPUT
      !        datjd     : date in julian days
      !
      !     OUTPUT
      !        obliquity : mean obliquity in degrees at datjd
      !******************************************************
      !
      !      use constantes
      implicit none

      real(kind = dp), intent(in) :: datjd
      real(kind = dp) :: obliquity

      real(kind = dp) :: t

      t = (datjd-(xjd2000)) / century ! parameters  referered to J2000

      !      obliquity =  84381.410d0/3600d0 -(46.80956*t + 0.000152*t*t+0.0019989*t*t*t)/3600d0
      !      obliquity =  84381.4118988d0/3600d0 -( 46.80972*t + 0.001998*t*t*t)/3600d0  ! used in the ephemeride version to LL for GAIA only (not the planets)

      obliquity = 84381.4119d0/3600d0 - (46.80956 * t + 0.000152 * t * t + 0.0019989 * t * t * t)/3600d0 !used for the GAIA Chebyshev

      return
   end function obliquity   
   !*****************************************************************
   subroutine dpsideps(xjd, dpsi, deps)
      !*****************************************************************
      !     Nutation in longitude and obliquity at the date t in JD
      !     dpsi and deps in arcsec
      !     Relatively mild accuracy  ~ 0.02 arcsec
      !     Use FM nutation_onemas(date1, date2, dpsi, deps) for 1 mas accuracy
	  !     F. Mignard. OCA/Cassiopee. March 2005
      !
      !      INPUT
      !         xjd : date in  julian days
      !
      !      OUTPUT
      !        dpsi : nutation in longitude in arcsec
      !        deps : nutation in obliquity in arcsec
      !*****************************************************************

      !     use constantes

      implicit none

      real(kind = dp), intent(in) :: xjd
      real(kind = dp), intent(out) :: dpsi, deps

      !     real(kind=dp)                 :: fp                  ! function fractional part
      real(kind = dp) :: t ! epoch from J2000 in centuries
      real(kind = dp) :: xml, xms, d, f, o ! Delaunay elements of the Moon


      !
      ! date since J2000 in julian centuries
      !
      t = (xjd-xjd2000) / century

      !
      !Moon mean anomaly at t  (round = 360d0 in module constantes)
      !
      xml = deupi * fp(0.374897d0 + 1325.552410d0 * t + 2.428d-5 * t * t + 3.98d-8 * t * t * t)
      !
      !Sun mean anomaly at t
      !
      xms = deupi * fp(0.993133d0 + 99.997361d0 * t - 4.26d-7 * t * t)
      !
      !Synodicl longitude of the Moon wrt the Sun
      !
      D = deupi * fp(0.827361d0 + 1236.853086d0 * t - 4.528d-6 * t * t)
      !
      !Latitude argument of the Moon
      !
      F = deupi * fp(0.259086d0 + 1342.227825d0 * t - 10.13d-6 * t * t - 7.9d-10 * t * t * t)

      !
      !Ascending node of the Moon
      !
      O = deupi * fp(0.347337d0 - 5.372605d0 * t)

      dpsi = (-17.1996d0 - 0.0174 * t) * sin(o) - 1.3187d0 * sin(2 * F - 2 * D + 2 * o) - 0.2274d0 * sin(2 * F + 2 * o) &
      +0.2062d0 * sin(2 * o) + 0.1426d0 * sin(xms) + 0.0712d0 * sin(xml) - 0.0517d0 * sin(xms + 2 * f - 2 * d + 2 * o) &
      -0.0386d0 * sin(2 * f + o) - 0.0301d0 * sin(xml + 2 * f + 2 * o) + 0.0217d0 * sin(-xms + 2 * f - 2 * d + 2 * o) &
      -0.0158d0 * sin(xml - 2 * d) + 0.0129d0 * sin(2 * f - 2 * d + o) + 0.0123d0 * sin(-xml + 2 * f + 2 * o)


      deps = 9.2025d0 * cos(o) + 0.5736d0 * cos(2 * F - 2 * D + 2 * o) + 0.0977d0 * cos(2 * F + 2 * o) - 0.0895d0 * cos(2 * o) &
      +0.0054d0 * cos(xms) - 0.0007 * cos(xml) + 0.0224 * cos(xms + 2 * f - 2 * d + 2 * o) + 0.0200d0 * cos(2 * f + o) &
      +0.0129d0 * cos(xml + 2 * f + 2 * o) - 0.0095d0 * cos(-xms + 2 * f - 2 * d + 2 * o) - 0.0001d0 * cos(xml - 2 * d) &
      -0.0070d0 * cos(2 * f - 2 * d + o) - 0.0053d0 * cos(-xml + 2 * f + 2 * o)

      return
 end subroutine dpsideps
   !******************************************************
   subroutine equtoecl(datjd, poseq, posecl)
      !******************************************************
      !
      !    Transformation of a vector from the equatorial frame to the ecliptic frame
      !
      !  Author  : F. Mignard OCA/CERGA
      !
      !
      !*** INPUT
      !       datjd   : JD of the date (used to compute the obliquity)
      !       poseq   : input vector in equatorial frame
      !
      !*** OUTPUT
      !       posecl  : output vector in ecliptic frame
      !
      !******************************************************
      !

      implicit none
      real(kind = dp), parameter :: phi = 0.05542d0 !Chapront et al. 2002 gamma to O_icrs in arcsec

      real(kind = dp), dimension(3), intent(in) :: poseq
      real(kind = dp), dimension(3), intent(out) :: posecl

      real(kind = dp) :: epsil, datjd, phid, posint(3)

      phid = phi/3600d0

      epsil = obliquity(datjd)

      call rota(poseq, -phid, 3, posint)

      call rota(posint, epsil, 1, posecl)

      return
   end subroutine equtoecl
   !******************************************************
   subroutine ecltoequ(datjd, posecl, poseq)
      !******************************************************
      !
      !    Transformation of a vector from the ecliptic frame to the equatorial frame
      !
      !  Author  : F. Mignard OCA/CERGA
      !
      !
      !*** INPUT
      !       datjd   : JD of the date (used to compute the obliquity)
      !       posecl  : input vector in ecliptic frame
      !
      !*** OUTPUT
      !       poseq   : output vector in equatorial frame
      !
      !******************************************************
      !

      implicit none
      real(kind = dp), parameter :: phi = 0.05542d0 !Chapront et al. 2002 gamma to O_icrs in arcsec
   !   real(kind=dp), parameter :: phi = 0.0d0  !Chapront et al. 2002 gamma to O_icrs in arcsec

      real(kind = dp), dimension(3), intent(in)  :: posecl
      real(kind = dp), dimension(3), intent(out) :: poseq

      real(kind = dp) :: epsil, datjd, phid, posint(3)

      phid = phi/3600d0
      epsil = obliquity(datjd)

      call rota(posecl, -epsil, 1, posint)

      call rota(posint, phid, 3, poseq) ! not include this statement for the current version of ephemeris for GAIA

      return
   end subroutine ecltoequ   
 
end module jmp_planetlib
