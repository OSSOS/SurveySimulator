module numutils

  use datadec
  use xvutils
  use rot

contains

  subroutine AppMag (r, delta, robs, h, g, alpha, mag, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes phase angle and apparent magnitude.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     r     : Sun-object distance (R8)
!     delta : Earth-object distance (R8)
!     robs  : Sun-Earth distance (R8)
!     h     : Absolute magnitude of object (R8)
!     g     : Slope of object (R8)
!
! OUTPUT
!     alpha : Phase angle (R8)
!     mag   : Apparent magnitude (R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : wrong input data
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) r
!f2py intent(in) delta
!f2py intent(in) robs
!f2py intent(in) h
!f2py intent(in) g
!f2py intent(out) alpha
!f2py intent(out) mag
!f2py intent(out) ierr
    implicit none

    integer :: ierr
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, raddeg = 180.0d0/Pi
    real (kind=8) :: r, delta, robs, h, g, alpha, mag, denom, phi1, phi2

    ierr = 0
    denom = 2.d0*r*delta
    if (denom .eq. 0.d0) then
       ierr = 10
       return
    end if
    alpha = dacos(dmin1((-robs**2 + delta**2 + r**2)/denom,1.d0))
    phi1 = exp(-3.33d0*(dtan(alpha/2.0d0))**0.63d0)
    phi2 = exp(-1.87d0*(dtan(alpha/2.0d0))**1.22d0)
    mag = 5.d0*dlog10(r*delta) + h - 2.5d0*dlog10((1.d0 - g)*phi1 + g*phi2)
!      write (6, '(7(f10.4,1x))')
!     $  r, delta, robs, alpha, h, mag,
!     $  2.5d0*dlog10((1.d0 - g)*phi1 + g*phi2)

    return
  end subroutine AppMag

  subroutine AbsMag (r, delta, robs, mag, g, alpha, h, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes phase angle and absolute magnitude.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : May 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     r     : Sun-object distance (R8)
!     delta : Earth-object distance (R8)
!     robs  : Sun-Earth distance (R8)
!     mag   : Apparent magnitude (R8)
!     g     : Slope of object (R8)
!
! OUTPUT
!     alpha : Phase angle (R8)
!     h     : Absolute magnitude of object (R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : wrong input data
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) r
!f2py intent(in) delta
!f2py intent(in) robs
!f2py intent(in) mag
!f2py intent(in) g
!f2py intent(out) alpha
!f2py intent(out) h
!f2py intent(out) ierr
    implicit none

    integer ierr
    real (kind=8) :: r, delta, robs, h, g, alpha, mag, mag0

    h = 0.d0
    call AppMag (r, delta, robs, h, g, alpha, mag0, ierr)
    if (ierr .ne. 0) then
       write (6, *) 'AppMag: something''s wrong !'
       write (6, *) 'ierr = :', ierr
       stop
    end if
    h = mag - mag0

    return
  end subroutine AbsMag

  subroutine dgauss (i, y)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives a random value with gaussian probability, with 0
! mean and standard deviation 1.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : January 1990
! Version 2 : January 2007
!             Modified to use RAN3 as random number generator rather
!             than PSALUN because it has a much longer periodicity.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     i     : Seed for random number generator (I4)
!
! OUTPUT
!     y     : Random value (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) i
!f2py intent(out) y
    implicit none

    integer :: i
    integer, save :: compte
    real (kind=8), save :: pi, x1, x2
    real (kind=8) :: y, y1, y2

    data &
         compte /0/, &
         pi /3.141592653589793238d0/

    if (compte.eq.0) then
       y1=ran3(i)
       y2=ran3(i)
       y1=dsqrt(-2.*dlog(y1))
       y2=2.*pi*y2
       x1=y1*dcos(y2)
       x2=y1*dsin(y2)
       compte=1
       y=x1
    else
       compte=0
       y=x2
    end if

    return
  end subroutine dgauss

  subroutine magran (mag_t, mag_er, seed, mag, magerr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine gives a randomized value of the magnitude given the
! theoretical magnitude and parameters to compute the uncertainty.
!
! Version 2
! This works for uncertainties given by the measurement on 1 frame only.
! Shouldn't try to combine several frame to estimate the error as this
! mostly account for zeropoint uncertainty and lightcurve.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : January 2006
! Version 2 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     mag_t : Theoretical magnitude of object (R8)
!     mag_er: Magnitude error parameters (6,n*R8)
!     seed  : Seed for random number generator (I4)
!
! OUTPUT
!     mag   : Randomized magnitude (R8)
!     magerr: Magnitude uncertainty (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) mag_t
!f2py intent(in) mag_er
!f2py intent(in,out) seed
!f2py intent(out) mag
!f2py intent(out) magerr

    implicit none

    real (kind=8) :: mag_t, mag_er(*), magerr, mag, tmp, mag_th
    integer :: seed, i

    mag_th = mag_t
!      tmp = log10(mag_er(2)/mag_er(1))/(mag_er(3)-21.d0)
    if (mag_th .le. 21.d0) then
       magerr = mag_er(1)
    else if (mag_th .le. mag_er(3)) then
!         magerr = mag_er(1)*10.d0**(tmp*(mag_th-21.d0))
       magerr = mag_er(1)*10.d0**(mag_er(2)*(mag_th-21.d0))
    else
!         magerr = max(mag_er(1)*10.d0**(tmp*(mag_er(3)-21.d0))
       magerr = max(mag_er(1)*10.d0**(mag_er(2)*(mag_er(3)-21.d0)) &
            - (mag_th - mag_er(3))*mag_er(4), 0.d0)
    end if
    call dgauss(seed, tmp)
    mag = mag_th + magerr*tmp
!      write (19, *) (mag_er(i), i=1,6)
!      write (19, *) mag_th, magerr, tmp, mag
    if (mag_th .gt. mag_er(5)) then
       mag = mag + (mag_th - mag_er(5))*mag_er(6)
    end if
!      write (19, *) mag
!      write (19, '(4(f6.3, 1x), i10)') mag_th, magerr, tmp, mag, seed

    return

  end subroutine magran

  subroutine LatLong (pos, long, lat, r)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine transforms cartesian coordinates into longitude,
! latitute and distance (almost spherical coordinates). If the input
! cordinates are in ICRF, then one obtains the RA and DEC
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2003
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     pos   : Object's cartesian coordinates
!
! OUTPUT
!     long  : Longitude (RA if ICRF)
!     lat   : Latitude (DEC if ICRF)
!     r     : Distance to center
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) pos
!f2py intent(out) long
!f2py intent(out) lat
!f2py intent(out) r
    implicit none

    type(t_v3d) :: pos
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi
    real (kind=8) :: lat, long, r

    r = dsqrt (pos%x**2 + pos%y**2 + pos%z**2)
    long = datan2(pos%y, pos%x)
    if (long .lt. 0.d0) long = long + TwoPi
    lat = asin(pos%z/r)

    return
  end subroutine LatLong

  subroutine RADECeclXV (pos, obspos, delta, ra, dec)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the RA and DEC of an object, defined by its
! barycentric ecliptic cartesian coordinates, with respect to an
! observatory, defined by its ICRF cartesian coordinates.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     pos   : Object barycentric ecliptic cartsian coordinates (3*R8)
!     obspos: Observatory ICRF cartsian coordinates (3*R8)
!
! OUTPUT
!     delta : Distance to observatory (R8)
!     ra    : Right Ascension (R8)
!     dec   : Declination (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) pos
!f2py intent(in) obspos
!f2py intent(out) delta
!f2py intent(out) ra
!f2py intent(out) dec
    implicit none

    type(t_v3d) :: pos, opos, obspos
    real (kind=8) :: ra, dec, delta
    integer :: ierr

! Compute ICRF cartesian coordinates
    call equat_ecl (-1, pos, opos, ierr)
    if (ierr .ne. 0) then
       write (6, *) 'Problem in conversion ecliptic -> equatorial'
    end if

! Compute RA and DEC
    opos%x = opos%x - obspos%x
    opos%y = opos%y - obspos%y
    opos%z = opos%z - obspos%z
    call LatLong (opos, ra, dec, delta)

    return
  end subroutine RADECeclXV

  FUNCTION ran3(idum)
!f2py intent(in,out) idum
    implicit none

    integer :: idum
    integer, parameter :: MBIG=1000000000, MSEED=161803398, MZ=0
    real (kind=8), parameter :: FAC=1.d0/MBIG
    integer, save :: iff, inext, inextp, ma(55)
    integer :: i, ii, mj, mk, k
    real (kind=8) :: ran3

    data iff /0/

    if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=abs(MSEED-abs(idum))
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
       end do
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
       end do
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran3=mj*FAC
    return
  end function ran3

  subroutine zero2pi (var)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This function resets variable 'var' to be between 0 and 2pi
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! B. Gladman  UBC
! Version 1 : January 2007
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT/OUPUT
!     var   : Variable to reset to be between 0 and 2*Pi (R8)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) var
    implicit none

! Calling arguments
    real (kind=8) :: var

!Some values better set up as parameters
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi

771 if (var .gt. TwoPi) then
       var = var - TwoPi
       goto 771
    endif
772 if (var .lt. 0.0d0) then
       var = var + TwoPi
       goto 772
    endif
    return
  end subroutine zero2pi

  subroutine cal2jul (iyyy, mm, dd, jul)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the Julian Day from time given in
! Year, Month, Day (decimal).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2003
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     iyyy  : Year (I4)
!     mm    : Month (I4)
!     dd    : Decimal Day (R8)
!
! OUTPUT
!     mjd   : Julian Day (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in) iyyy
!f2py intent(in) mm
!f2py intent(in) dd
!f2py intent(out) jul
    implicit none

    integer :: mm, iyyy, id, juld
    real (kind=8) :: dd, idfrac, jul

    id = INT(dd)
    idfrac = dd - id
    juld = julday(mm, id, iyyy)
    if (idfrac .lt. 0.5d0) then
       juld = juld - 1
       idfrac = idfrac + 0.5d0
    else
       idfrac = idfrac - 0.5d0
    end if
    jul = juld + idfrac

    return
  end subroutine cal2jul

  function julday(MM,ID,IYYY)

    implicit none

    integer, parameter :: igreg=15+31*(10+12*1582)
    integer :: mm, id, iyyy, jy, jm, ja, julday

    IF (IYYY.EQ.0) stop 'There is no Year Zero.'
    IF (IYYY.LT.0) IYYY=IYYY+1
    IF (MM.GT.2) THEN
       JY=IYYY
       JM=MM+1
    ELSE
       JY=IYYY-1
       JM=MM+13
    ENDIF
    JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
    IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
       JA=INT(0.01*JY)
       JULDAY=JULDAY+2-JA+INT(0.25*JA)
    ENDIF
    RETURN
  END function julday

  subroutine ObjAbs (o_p, jday, mag, code, gb, alpha, h, ra, dec, ierr)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine determines the absolute magnitude of an object and its
! sky position (RA, DEC) given its orbital elements (Berstein &
! Kushalani format), measured magnitude and epoch of observation.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : May 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     o_p   : Orbital elements with time of peri (t_orb_p)
!     jday  : Epoch of observation [JD] (R8)
!     mag   : Apparent magnitude of object (R8)
!     code  : Observatory code (I4)
!              001 : GAIA
!              002 : Geocentric, Mignard's code
!              500 : Geocentric
!     gb    : opposition surge factor, Bowell formalism (R8)
!
! OUTPUT
!     alpha : Phase angle [rad] (R8)
!     h     : Absolute magnitude of object (R8)
!     ra    : Right Ascension (R8)
!     dec   : Declination (R8)
!     ierr  : Error code (I4)
!                0 : nominal run
!               10 : wrong input data
!              100 : date of call earlier than xjdbeg
!              200 : date of call later   than xjdend
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) o_p
!f2py intent(in) jday
!f2py intent(in) mag
!f2py intent(in) gb
!f2py intent(out) alpha
!f2py intent(out) h
!f2py intent(out) ra
!f2py intent(out) dec
!f2py intent(out) ierr
    implicit none

    type(t_orb_p) :: o_p
    type(t_orb_m) :: o_m
    type(t_v3d) :: pos, obpos, tmp
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, &
         drad = Pi/180.0d0
    integer, parameter :: screen = 6
    real (kind=8) :: h, jday, alpha, ra, dec, gb, mag
    real (kind=8), save :: r, delta, ros
    integer :: ierr, code

    call ObsPos (code, jday, obpos, tmp, ros, ierr)
    if (ierr .ne. 0) then
       write (screen, *) 'Error while computing observatory''s position.'
       write (screen, *) 'ierr = ', ierr
       return
    end if
    o_m%a = o_p%a
    o_m%e = o_p%e
    o_m%inc = o_p%inc
    o_m%node = o_p%node
    o_m%peri = o_p%peri
    o_m%m = (twopi/(o_p%a**1.5d0*365.25d0))*(jday-o_p%tperi)
    o_m%m = o_m%m - int(o_m%m/twopi)*twopi
    call pos_cart (o_m, pos)
    call DistSunEcl (jday, pos, r)
    call RADECeclXV (pos, obpos, delta, ra, dec)
    call AbsMag (r, delta, ros, mag, gb, alpha, h, ierr)
    if (ierr .ne. 0) then
       write (screen, *) 'AbsMag: something''s wrong !'
       write (screen, *) 'ierr = :', ierr
       return
    end if

    return
  end subroutine ObjAbs

end module numutils
