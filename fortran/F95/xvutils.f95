module xvutils

  use datadec
  use DE_utils
  use DE_const
  use DE_ephhdr
  use rot
  use ioutils

contains

  subroutine DistSunEcl (jday, pos, r)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes distance from object defined by barycentric
! ecliptic coordinates to Sun at given jday.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : February 2023
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
    type(t_v3d), save :: pos_s, vel_s
    integer :: ierr
    real (kind=8), save :: jday_s
    real (kind=8) :: pv(6)
    logical, save :: done

    data done /.false./

    if ((.not. done) .or. (jday .ne. jday_s)) then
! Get (Equatorial) barycentric position and velocity of Sun
       call pleph(jday, idtarg(1), nctr, pv)
       if (pv(1) .le. -9d99) stop 'Failed while getting barycentric position of Sun.'
       pos_s%x = pv(1)
       pos_s%y = pv(2)
       pos_s%z = pv(3)
       vel_s%x = pv(4)
       vel_s%y = pv(5)
       vel_s%z = pv(6)
! Compute Ecliptic cartesian coordinates
       call equat_ecl (1, pos_s, pos_s, ierr)
       if (ierr .ne. 0) then
          write (6, *) 'Problem in conversion equatorial -> ecliptic'
       end if
       done = .true.
       jday_s = jday
    end if

    r = dsqrt((pos%x - pos_s%x)**2 + (pos%y - pos_s%y)**2 &
         + (pos%z - pos_s%z)**2)
    return

  end subroutine DistSunEcl

  subroutine topo (geolong, geolat, height, pos)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the geocentric coordinates from the geodetic
! (standard map-type) longitude, latitude, and height. These are
! assumed to be in radians, radians and meters respectively.
! Notation generally follows 1992 Astr Almanac, p. K11
! Returns the cartesian coordinates of the observatory.
! Units : AU
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     geolong: longitude of observatory (rd) (R8)
!     geolat: latitude of observatory (rd) (R8)
!     height: altitude of observatory above ground (m) (R8)
!
! OUTPUT
!     pos   : Cartesian coordinates of observatory (AU) (3*R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) geolong
!f2py intent(in) geolat
!f2py intent(in) height
!f2py intent(out) pos
    implicit none

    real (kind=8), intent(in) :: geolong, geolat, height
    real (kind=8), dimension(3), intent(out) :: pos
    real (kind=8) :: denom, C_geo, S_geo

    denom = (1.0d0 - Flatten) * sin(geolat)
    denom = cos(geolat)**2 + denom*denom
    C_geo = 1.0d0 / sqrt(denom)
    S_geo = (1.0d0 - Flatten)**2 * C_geo
    C_geo = C_geo* Equat_Rad + height ! deviation from almanac notation -- include height here.
    S_geo = S_geo* Equat_Rad + height
    pos(1) = C_geo * cos(geolat) * cos(geolong)
    pos(2) = C_geo * cos(geolat) * sin(geolong)
    pos(3) = S_geo * sin(geolat);

! convert to AU, keeping in mind that Horizons was km and this is m
    pos(1) = pos(1)/(1000.0d0*km2AU)
    pos(2) = pos(2)/(1000.0d0*km2AU)
    pos(3) = pos(3)/(1000.0d0*km2AU)

    return
  end subroutine topo

  subroutine observatory_geocenter (code, t, pos, obfile)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the cartesian coordinates of the observatory
! with respect to the geocenter, at time t.
! Returns cartisian coordiantes.
! Reference frame : ICRF
! Units : AU
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     code  : Observatory code (I4)
!              000 : Solar System Barycenter
!              001 : GAIA
!              002 : Geocentric, Mignard's code
!              500 : Geocentric
!              <0  : -lun of a file handle to read JPL state vectors from
!     t     : Time of observation (Julian day, not MJD) (R8)
!     obfile: Name of file with observatory list (optional) (CH)
!
! OUTPUT
!     pos   : Cartesian coordinates of observatory (AU) (3*R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) code
!f2py intent(in) t
!f2py intent(out) pos
    implicit none

    integer (kind=4), intent(in) :: code
    real (kind=8), intent(in) :: t
    character (len=*), intent(in), optional :: obfile
    real (kind=8), dimension(3), intent(out) :: pos
    integer (kind=4) :: i
    real (kind=8) :: obslon, obslat, obsalt, obslmst

    if (nsites .lt. 0) then
       call read_observatories(obfile)
    end if

    i = 1
    do while ((i .le. nsites) .and. (code .ne. sitelist(i)%code))
       i = i + 1
    end do
    if (i .gt. nsites) then
       print *, 'ERROR, unknown observatory code', code
       stop 1
    end if

    obslon = sitelist(i)%lon
    obslat = sitelist(i)%lat
    obsalt = sitelist(i)%altitude

    obslmst = lst(t, obslon)
    call topo (obslmst/12.0d0*Pi, obslat, obsalt, pos)

    return
  end subroutine observatory_geocenter

  subroutine ObsPos (code, t, pos, vel, r, ierr, ephfile, obfile)
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
!     ephfile: Ephemerides file name (optional) (CH)
!     obfile: Name of file with observatory list (optional) (CH)
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
    character (len=*), intent(in), optional :: ephfile
    character (len=*), intent(in), optional :: obfile
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: t
    real (kind=8), intent(out) :: r
    type(t_v3d) :: pos_s, vel_s
    integer (kind=4) :: istat, ntar, i1, i2
    data ntar /4/ ! Earth index in idtarg
    real (kind=8) :: v_pos(3), pv(6)
    logical :: first_obspos, finished
    data first_obspos /.true./
    save first_obspos, ntar

    if (first_obspos) then
! Set epehmerides file name, if need be
       if (present(ephfile)) then
! Initialize ephemerides
          call init_DE(ephfile)
       else
! Initialize ephemerides
          call init_DE()
       end if
       first_obspos = .false.
    end if

    ierr = 0
! Get (Equatorial) barycentric position and velocity of Sun
    call pleph(t, idtarg(1), nctr, pv)
    if (pv(1) .le. -9d99) stop 'Failed while getting barycentric position of Sun.'
    pos_s%x = pv(1)
    pos_s%y = pv(2)
    pos_s%z = pv(3)
    vel_s%x = pv(4)
    vel_s%y = pv(5)
    vel_s%z = pv(6)
    if (code .eq. 0) then
       pos%x = 0.0d0
       pos%y = 0.0d0
       pos%z = 0.0d0
       vel%x = 0.0d0
       vel%y = 0.0d0
       vel%z = 0.0d0
    else if (code .eq. 1) then
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
        
! Get barycentric equatorial Earth coordinates (ICRF).
       call pleph(t, idtarg(ntar), nctr, pv)
       if (pv(1) .le. -9d99) stop 'Failed while getting barycentric position of Earth.'
       pos%x = pv(1)
       pos%y = pv(2)
       pos%z = pv(3)
       vel%x = pv(4)
       vel%y = pv(5)
       vel%z = pv(6)
       if (code .ne. 500) then
          call observatory_geocenter(code, t, v_pos, obfile)
          pos%x = pos%x + v_pos(1)
          pos%y = pos%y + v_pos(2)
          pos%z = pos%z + v_pos(3)
       end if
    end if
! Finally, computes distance from observatory to Sun.
    r = dsqrt((pos%x - pos_s%x)**2 + (pos%y - pos_s%y)**2 &
         + (pos%z - pos_s%z)**2)

  end subroutine ObsPos

end module xvutils
