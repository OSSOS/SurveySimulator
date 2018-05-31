module gimeobjut

  use datadec
  use numutils
  use ioutils

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Generic model routines for Survey Simulator, version 2.0 for OSSOS
!
! Calling sequence in SurveySimulator.f (survey simulator driver) is:
!
! Loop (some condition):
!     call GiMeObj(arg_list_1)
!     Check model ended:
!         set exit condition
!     call Detos1(arg_list_2)
!     Check detection and tracking:
!         store results
!
! where arg_list_1 is
! (filena, seed, a, e, inc, node, peri, M, epoch, h, color, gb, ph,
!  period, amp, comp, ierr)
! with:
!
! INPUT
!     filena: name of containing description of model, read in by the
!             model subroutine "modname" (CH)
!     seed  : Random number generator seed (I4)
!
! OUTPUT
!     a     : semimajor axis (R8)
!     e     : eccentricity (R8)
!     inc   : Inclination [rad] (R8)
!     node  : Longitude of node [rad] (R8)
!     peri  : Argument of perihelion [rad] (R8)
!     M     : Mean anomaly [rad] (R8)
!     epoch : epoch of the orbital elements, in Julian Day (R8)
!     h     : absolute magnitude of object in band filter "x" (R8)
!     color : array of colors "y-x", where the index of "y" is as
!             described in detos1 (10*R8)
!                color(1) : g-x
!                color(2) : r-x
!                color(3) : i-x
!                color(4) : z-x
!                color(5) : u-x
!                color(6) : V-x
!                color(7) : B-x
!                color(8) : R-x
!                color(9) : I-x
!     gb    : opposition surge factor, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The model subroutines can access files using logical unit numbers from
! 20 upward. This range in reseved for them and won't be used by the
! drivers nor SurveySubs routines.
!
! It is good practice that when first started, the GiMeObj routine
! writes a file describing the model used, the versino and the date of
! the routine.
!
! Since this routine is called once for every object created, it needs
! to get all the required parameters once when it is called the first
! time, then save these values for future use.
!
! The following routine gives a working example of a model routine. It
! is probably worth reading it through.
!
! The survey simulator expects orbital elements with respect to ecliptic
! reference frame.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File generated on 2013-07-01
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (filena, seed, o_m, epoch, h, color, gb, ph, period, &
       amp, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model stored in a file.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : May 2013
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     filena: name of containing description of model, read in by the
!             model subroutine "modname" (CH)
!     seed  : Random number generator seed (I4)
!
! OUTPUT
!     o_m   : orbital elements of object (t_orb_m)
!     epoch : Time of elements [JD] (R8)
!     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
!     color : Array of colors (10*R8)
!                colors(1) : g-x
!                colors(2) : r-x
!                colors(3) : i-x
!                colors(4) : z-x
!                colors(5) : u-x
!                colors(6) : V-x
!                colors(7) : B-x
!                colors(8) : R-x
!                colors(9) : I-x
!     gb    : opposition surge factor, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code (I4)
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The user can fill the 100-character 'commen' string any way they
! wish; this comment string will be printed in the driver on the output
! line of each detection.  Examples of the comment might be resonance name
! and libration amplitude, or the name of a component in the GiMeObj model
! that the object responds to.  The nchar variable (passed back to Driver)
! allows the user to pring only the first nchar characters of this string.
!
! This routine uses logical unit 20 to access the file containing the model.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) filena
!f2py intent(in,out) seed
!f2py intent(out) o_m
!f2py intent(out) epoch
!f2py intent(out) h
!f2py intent(out) color
!f2py intent(out) gb
!f2py intent(out) ph
!f2py intent(out) period
!f2py intent(out) amp
!f2py intent(out) commen
!f2py intent(out) nchar
!f2py intent(out) ierr
    implicit none

! Calling arguments
    integer :: ierr, seed, nchar
    type(t_orb_m) :: o_m
    real (kind=8) :: epoch, h, color(10), gb, ph, period, amp
    character :: filena*(*), commen*100

! Some values better set up as parameters
    integer, parameter :: &
         n_obj_max = 100,       &! Maximum number of objects we can read at once
         lun_d = 20,            &! Logical unit number for data file reading
         lun_ll = 21,           &! Logical unit number for logging
         screen = 6              ! Logical unit number for screen writing
    real (kind=8), parameter :: &
         Pi = 3.141592653589793d0, &! Pi
         TwoPi = 2.0d0*Pi,         &! 2*Pi
         drad = Pi/180.0d0          ! Degree to radian convertion: Pi/180

! Internal storage
    type(t_orb_m), save :: &
         obj_o(n_obj_max)        ! Array of orbital elements
    real (kind=8), save :: &
         obj_h(n_obj_max),      &! Array of absolute magnitudes of object
                                 ! in 'x' band
         obj_jday(n_obj_max),   &! Array of times of elements [JD]
         random,                &! Random number
         color0(10)              ! Color parameters of model
    character, save :: &
         comp(n_obj_max)*10      ! Array of strings telling the component
                                 ! the object belongs too
    integer, save :: &
         i_obj,                 &! Index of current object in array
         n_obj,                 &! Number of objects in arrays
         ierr_d,                &! Return code for GetDistrib
         i1, i2                  ! Dummy indices
    logical, save :: &
         end_of_file,           &! Tells if we've reached end of file
         finished,              &! Tells if string is empty
         first                   ! Tells if first call to routine

! Lightcurve and opposition surge effect parameters
    real (kind=8), save :: gb0, ph0, period0, amp0

! Sets initial values
    data &
         i_obj / 0/,            &! Starts from object 0
         end_of_file /.false./, &! We are not at the end yet
         first /.true./,        &! First call
         gb0     /-0.12d0/,     &! Opposition surge effect
         ph0     / 0.00d0/,     &! Initial phase of lightcurve
         period0 / 0.60d0/,     &! Period of lightcurve
         amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

! This is the first call
    if (first) then
       call read_file_name(filena, i1, i2, finished, len(filena))
! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a,a,a)') &
            '# Model from file ', filena(i1:i2), ', version 1.0, 2013-05-14'
       close (lun_ll)
! Change "first" so this is not called anymore
       first = .false.
    end if

! Check if there are still objects available in the arrays
    if (i_obj .le. 0) then
! Get new objects. If this is the first call to GetDistrib, it will
! first open the data file, otherwise, will simply return the following
! objects. When reaching the end of the file, returns "ierr_d = 30".
       call GetDistrib (filena, lun_d, n_obj_max, n_obj, obj_o, &
            obj_h, obj_jday, color0, comp, ierr_d)
       if (n_obj .le. 0) then
          ierr = -20
          return
       end if
       if (ierr_d .ne. 0) then
          if (ierr_d .eq. 10) then
             write (screen, *) 'Unable to open ', filena
          else if (ierr_d .eq. 30) then
             end_of_file = .true.
             goto 100
          else
             write (screen, *) 'Unknown return code in read_obj:', ierr_d
          end if
! If we get here, there is something really wrong, better return with
! panic code.
          ierr = -20
          return
       end if
    end if
! Ok, now we have data to send back
100 continue

! Now extract an object from the arrays
    i_obj = i_obj + 1

! Determine the object elements (so they're the same for all surveys).
    o_m = obj_o(i_obj)

! If the angles are "not defined" (i.e. < -twopi) they are drawn at
! random. If jday is "not defined" (i.e. < 0.), mean anomaly is drawn
! at random.
    if (o_m%node .lt. -twopi) then
       random = ran3(seed)
       o_m%node = random*twopi
    end if
    if (o_m%peri .lt. -twopi) then
       random = ran3(seed)
       o_m%peri = random*twopi
    end if
    epoch = obj_jday(i_obj)
    if ((o_m%m .lt. -twopi) .or. (epoch .lt. 0.d0)) then
       random = ran3(seed)
       o_m%m = random*twopi
       epoch = 2453157.5d0
    end if

! Get the absolute magnitude H.
    h = obj_h(i_obj)

! Copy the component name
    commen = comp(i_obj)
    nchar = len(comp(i_obj))

! Define values for lightcurve and opposition surge effect
    gb = gb0
    ph = ph0
    period = period0
    amp = amp0

! Get colors for object
    do ierr_d = 1, 10
       color(ierr_d) = color0(ierr_d)
    end do

! Prepare return code
    ierr = 0
    if (i_obj .eq. n_obj) then
       i_obj = 0
       if (end_of_file) then
! We've reached the end of the file, better tell the caller.
          ierr = 100
       end if
    end if

    return
  end subroutine GiMeObj

  subroutine GetDistrib (distri, lun_d, n_max, n_obj, obj_o, obj_h, obj_t, &
       color, comp, ierr)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine reads in an orbit distribution file.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     distri: Distribution file name (CH)
!     lun_d : Logical unit for file (I4)
!     n_max : Maximum number of objects to read (I4)
!
! OUTPUT
!     n_obj : Number of objects read (I4)
!     obj_o : Arry of orbital elements (n*t_orb_m)
!     obj_h : Absolute magnitude of object (n*R8)
!     obj_t : Time of elements of object (n*R8)
!     color : Array of colors (10*R8)
!     comp  : Dynamical componant the object belongs to (n*CH)
!     g     : Slope of object (R8)
!     ierr  : Error code (I4)
!                0 : nominal run, reached maximum number of objects
!               10 : unable to open filen
!               20 : error reading record
!               30 : end of file reached
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    implicit none

    integer :: n_max
    type(t_orb_m) :: obj_o(:), o_m
    real (kind=8) :: obj_h(:), obj_t(:), color(:), h, g, mag, r, alpha
    integer :: lun_d, n_obj, ierr
    character :: distri*(*), comp(*)*10, co*10
    real (kind=8), save :: jday

    ierr = 0

! Hard coded slope for magnitude ! Bad boy !
    g = -0.12d0

! Open and read in object distribution
    n_obj = 0
100 continue
    call read_obj (distri, lun_d, o_m, h, jday, color, co, ierr)

    if (ierr .ne. 0) then
       if (ierr .eq. 10) then
          write (6, *) 'Unable to open ', distri
       else if (ierr .eq. 20) then
          write (6, *) 'Error reading ', distri
          write (6, *) 'Object number: ', n_obj
          goto 100
       else if (ierr .eq. 30) then
          goto 110
       else
          write (6, *) 'Unknown return code in read_obj.'
       end if
       return
    end if

    n_obj = n_obj + 1
    obj_o(n_obj) = o_m
    obj_h(n_obj) = h
    obj_t(n_obj) = jday
    if (co .ne. '          ') comp(n_obj) = co
    if (n_obj .ge. n_max) then
       return
    end if
    goto 100
110 continue

    return
  end subroutine GetDistrib

  subroutine read_obj (filen, lun_in, o_m, h, jday, color, co, ierr)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine opens and reads in the object element file.
! Angles are returned in radian.
! Potentially use a common to return the H distribution parameters.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : For L7 data release, June 2010
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     filen : object element file name
!     lun_in: File unit
!
! OUTPUT
!     o_m   : orbital elements (t_orb_m)
!     h     : Absolute magnitude (R8)
!     jday  : Time of elements (R8)
!     color : Array of colors (10*R8)
!     co    : Dynamical componant the object belongs to (CH)
!     ierr  : Error code
!                0 : nominal run
!               10 : unable to open filen
!               20 : error reading record
!               30 : end of file reached
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    implicit none

    type(t_orb_m) :: o_m
    real (kind=8) :: h, jday, color(*)
    real (kind=8), save :: jd
    real (kind=8), parameter :: Pi = 3.141592653589793d0, drad = Pi/180.0D0
    integer, parameter :: nw_max = 20
    integer :: lun_in, ierr, j, nw, lw(nw_max)
    character :: line*100, filen*(*), word(nw_max)*80, co*10
    logical, save :: opened

    data opened /.false./

    ierr = 0
    if (.not. opened) then
       open (unit=lun_in, file=filen, status='old', err=1000)
       opened = .true.
       jd = -1.d0
       color(1:10) = 0.d0
    end if

1500 continue
    do j = 1, len(line)
       line(j:j) = ' '
    end do
    read (lun_in, '(a)', err=2000, end=3000) line
    if (line(1:1) .eq. '#') then
       if (line(1:25) .eq. '# Epoch of elements: JD =') then
          read (line(26:100), *, err=1500, end=1500) jd
       end if
       if (line(1:10) .eq. '# Colors =') then
          read (line(11:), *, err=1500, end=1500) (color(j),j=1,10)
       end if
       goto 1500
    end if
    jday = jd
    call parse (line, nw_max, nw, word, lw)
    if (nw .lt. 7) goto 2000
    read (word(1), *) o_m%a
    read (word(2), *) o_m%e
    read (word(3), *) o_m%inc
    o_m%inc = o_m%inc*drad
    read (word(4), *) o_m%node
    read (word(5), *) o_m%peri
    read (word(6), *) o_m%m
    read (word(7), *) h
    o_m%node = o_m%node*drad
    o_m%peri = o_m%peri*drad
    o_m%m = o_m%m*drad
    co = '          '
    if (nw .ge. 9) co = word(9)
    return

1000 continue
    ierr = 10
    return

2000 continue
    ierr = 20
    return

3000 continue
    ierr = 30
    close (lun_in)
    opened = .false.
    return

  end subroutine read_obj

end module gimeobjut
