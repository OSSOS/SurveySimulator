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
    type(t_v3d) :: pos_s, vel_s
    integer (kind=4) :: istat, ntar
    data ntar /4/ ! Earth index in idtarg
    real (kind=8) :: v_pos(3), pv(6)
    logical :: first_obspos
    data first_obspos /.true./
    save first_obspos, ntar

    if (first_obspos) then
! Initialize ephemerides
       call init_DE()
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
          call observatory_geocenter(code, t, v_pos)
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
