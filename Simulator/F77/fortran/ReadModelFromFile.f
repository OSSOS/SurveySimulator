c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Generic model routines for Survey Simulator, version 2.0 for OSSOS
c
c Calling sequence in SurveySimulator.f (survey simulator driver) is:
c
c Loop (some condition):
c     call GiMeObj(arg_list_1)
c     Check model ended:
c         set exit condition
c     call Detos1(arg_list_2)
c     Check detection and tracking:
c         store results
c
c where arg_list_1 is
c (filena, seed, a, e, inc, node, peri, M, epoch, h, color, gb, ph,
c  period, amp, comp, ierr)
c with:
c
c INPUT
c     filena: name of containing description of model, read in by the
c             model subroutine "modname" (CH)
c     seed  : Random number generator seed (I4)
c
c OUTPUT
c     a     : semimajor axis (R8)
c     e     : eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     M     : Mean anomaly [rad] (R8)
c     epoch : epoch of the orbital elements, in Julian Day (R8)
c     h     : absolute magnitude of object in band filter "x" (R8)
c     color : array of colors "y-x", where the index of "y" is as
c             described in detos1 (10*R8)
c                color(1) : g-x
c                color(2) : r-x
c                color(3) : i-x
c                color(4) : z-x
c                color(5) : u-x
c                color(6) : V-x
c                color(7) : B-x
c                color(8) : R-x
c                color(9) : I-x
c     gb    : opposition surge factor, Bowell formalism (R8)
c     ph    : phase of lightcurve at epoch [rad] (R8)
c     period: period of lightcurve [day] (R8)
c     amp   : amplitude of lightcurve [mag] (R8)
c     commen: user specified string containing whatever the user wants (CH*100)
c     nchar : number of characters in the comment string that should be
c             printed out in output files if the object is detected;
c             maximum of 100 (I4)
c     ierr  : return code
c                  0 : nominal run, things are good
c                100 : end of model, exit after checking this object
c                -10 : could not get all orbital elements, skip object
c                -20 : something went grossly wrong, should quit
c
c The model subroutines can access files using logical unit numbers from
c 20 upward. This range in reseved for them and won't be used by the
c drivers nor SurveySubs routines.
c
c It is good practice that when first started, the GiMeObj routine
c writes a file describing the model used, the versino and the date of
c the routine.
c
c Since this routine is called once for every object created, it needs
c to get all the required parameters once when it is called the first
c time, then save these values for future use.
c
c The following routine gives a working example of a model routine. It
c is probably worth reading it through.
c
c The survey simulator expects orbital elements with respect to ecliptic
c reference frame.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File generated on 2013-07-01
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      subroutine GiMeObj (filena, seed, a, e, inc, node, peri, M,
     $  epoch, h, color, gb, ph, period, amp, commen, nchar, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c This routine generates an object from a model stored in a file.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2013
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT
c     filena: name of containing description of model, read in by the
c             model subroutine "modname" (CH)
c     seed  : Random number generator seed (I4)
c
c OUTPUT
c     a     : semimajor axis (R8)
c     e     : eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     M     : Mean anomaly [rad] (R8)
c     epoch : Time of elements [JD] (R8)
c     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
c     color : Array of colors (10*R8)
c                colors(1) : g-x
c                colors(2) : r-x
c                colors(3) : i-x
c                colors(4) : z-x
c                colors(5) : u-x
c                colors(6) : V-x
c                colors(7) : B-x
c                colors(8) : R-x
c                colors(9) : I-x
c     gb    : opposition surge factor, Bowell formalism (R8)
c     ph    : phase of lightcurve at epoch [rad] (R8)
c     period: period of lightcurve [day] (R8)
c     amp   : amplitude of lightcurve [mag] (R8)
c     commen: user specified string containing whatever the user wants (CH*100)
c     nchar : number of characters in the comment string that should be
c             printed out in output files if the object is detected;
c             maximum of 100 (I4)
c     ierr  : return code (I4)
c                  0 : nominal run, things are good
c                100 : end of model, exit after checking this object
c                -10 : could not get all orbital elements, skip object
c                -20 : something went grossly wrong, should quit
c
c The user can fill the 100-character 'commen' string any way they
c wish; this comment string will be printed in the driver on the output
c line of each detection.  Examples of the comment might be resonance name
c and libration amplitude, or the name of a component in the GiMeObj model
c that the object responds to.  The nchar variable (passed back to Driver)
c allows the user to pring only the first nchar characters of this string.
c
c This routine uses logical unit 20 to access the file containing the model.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) filena
Cf2py intent(in,out) seed
Cf2py intent(out) a
Cf2py intent(out) e
Cf2py intent(out) inc
Cf2py intent(out) node
Cf2py intent(out) peri
Cf2py intent(out) M
Cf2py intent(out) epoch
Cf2py intent(out) h
Cf2py intent(out) color
Cf2py intent(out) gb
Cf2py intent(out) ph
Cf2py intent(out) period
Cf2py intent(out) amp
Cf2py intent(out) commen
Cf2py intent(out) nchar
Cf2py intent(out) ierr

      implicit none

c Calling arguments
      integer*4 ierr, seed, nchar
      real*8 a, e, inc, node, peri, M, epoch, h, color(10), gb, ph,
     $  period, amp
      character filena*(*), commen*100

c Some values better set up as parameters
      integer*4
     $  n_obj_max,              ! Maximum number of objects we can read at once
     $  lun_d,                  ! Logical unit number for data file reading
     $  lun_ll,                 ! Logical unit number for logging
     $  screen                  ! Logical unit number for screen writing
      real*8
     $  Pi,                     ! Pi
     $  TwoPi,                  ! 2*Pi
     $  drad                    ! Degree to radian convertion: Pi/180
c Set the values
      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, lun_d = 20,
     $  drad = Pi/180.0d0, n_obj_max = 1000, screen = 6, lun_ll = 21)

c Internal storage
      real*8
     $  obj_a(n_obj_max),       ! Array of semimajor axis
     $  obj_e(n_obj_max),       ! Array of eccentricities
     $  obj_i(n_obj_max),       ! Array of inclinations [rad]
     $  obj_node(n_obj_max),    ! Longitude of nodes [rad]
     $  obj_peri(n_obj_max),    ! Array of arguments of perihelion [rad]
     $  obj_m(n_obj_max),       ! Array of mean anomalies [rad]
     $  obj_h(n_obj_max),       ! Array of absolute magnitudes of object
                                ! in 'x' band
     $  obj_jday(n_obj_max),    ! Array of times of elements [JD]
     $  random,                 ! Random number
     $  ran3,                   ! Random number generator
     $  color0(10)              ! Color parameters of model
      character
     $  date*8,                 ! Date of execution
     $  time*10,                ! Time of execution
     $  zone*5,                 ! Time zone
     $  comp(n_obj_max)*10      ! Array of strings telling the component
                                ! the object belongs too
      integer*4
     $  values(8),              ! Values of date and time
     $  i_obj,                  ! Index of current object in array
     $  n_obj,                  ! Number of objects in arrays
     $  ierr_d,                 ! Return code for GetDistrib
     $  i1, i2                  ! Dummy indices
      logical
     $  end_of_file,            ! Tells if we've reached end of file
     $  finished,               ! Tells if string is empty
     $  first                   ! Tells if first call to routine

c Lightcurve and opposition surge effect parameters
      real*8
     $  gb0, ph0, period0, amp0

c Ran3 is a function defined outside this routine
      external ran3

c Sets initial values
      data
     $  i_obj / 0/,             ! Starts from object 0
     $  end_of_file /.false./,  ! We are not at the end yet
     $  first /.true./,         ! First call
     $  gb0     /-0.12d0/,      ! Opposition surge effect
     $  ph0     / 0.00d0/,      ! Initial phase of lightcurve
     $  period0 / 0.60d0/,      ! Period of lightcurve
     $  amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

c Make sure variables are retained from call to call
      save

c This is the first call
      if (first) then
         call read_file_name(filena, i1, i2, finished, len(filena))
c Writes a file describing the model that was used.
         open (unit=lun_ll, file='ModelUsed.dat', access='sequential',
     $     status='unknown')
         write (lun_ll, '(a,a,a)')
     $     '# Model from file ', filena(i1:i2),
     $     ', version 1.0, 2013-05-14'
         call date_and_time(date, time, zone, values)
         write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ',
     $     date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' 
     $     //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
         close (lun_ll)
c Change "first" so this is not called anymore
         first = .false.
      end if

c Check if there are still objects available in the arrays
      if (i_obj .le. 0) then
c Get new objects. If this is the first call to GetDistrib, it will
c first open the data file, otherwise, will simply return the following
c objects. When reaching the end of the file, returns "ierr_d = 30".
         call GetDistrib (filena, lun_d, n_obj_max, n_obj, obj_a,
     $     obj_e, obj_i, obj_node, obj_peri, obj_m, obj_h, obj_jday,
     $     color0, comp, ierr_d)
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
               write (screen, *) 'Unknown return code in read_obj:',
     $           ierr_d
            end if
c If we get here, there is something really wrong, better return with
c panic code.
            ierr = -20
            return
         end if
      end if
c Ok, now we have data to send back
 100  continue

c Now extract an object from the arrays
      i_obj = i_obj + 1

c Determine the object elements (so they're the same for all surveys).
      a = obj_a(i_obj)
      e = obj_e(i_obj)
      inc = obj_i(i_obj)

c If the angles are "not defined" (i.e. < -twopi) they are drawn at
c random. If jday is "not defined" (i.e. < 0.), mean anomaly is drawn
c at random.
      node = obj_node(i_obj)
      if (node .lt. -twopi) then
         random = ran3(seed)
         node = random*twopi
      end if
      peri = obj_peri(i_obj)
      if (peri .lt. -twopi) then
         random = ran3(seed)
         peri = random*twopi
      end if
      M = obj_m(i_obj)
      epoch = obj_jday(i_obj)
      if ((M .lt. -twopi) .or. (epoch .lt. 0.d0)) then
         random = ran3(seed)
         M = random*twopi
         epoch = 2453157.5d0
      end if

c Get the absolute magnitude H.
      h = obj_h(i_obj)

c Copy the component name
      commen = comp(i_obj)
      nchar = len(comp(i_obj))

c Define values for lightcurve and opposition surge effect
      gb = gb0
      ph = ph0
      period = period0
      amp = amp0

c Get colors for object
      do ierr_d = 1, 10
         color(ierr_d) = color0(ierr_d)
      end do

c Prepare return code
      ierr = 0
      if (i_obj .eq. n_obj) then
         i_obj = 0
         if (end_of_file) then
c We've reached the end of the file, better tell the caller.
            ierr = 100
         end if
      end if

      return
      end

      subroutine GetDistrib (distri, lun_d, n_max, n_obj, obj_a, obj_e,
     $  obj_i, obj_no, obj_pe, obj_m, obj_h, obj_t, color, comp,
     $  ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine reads in an orbit distribution file.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     distri: Distribution file name (CH)
c     lun_d : Logical unit for file (I4)
c     n_max : Maximum number of objects to read (I4)
c
c OUTPUT
c     n_obj : Number of objects read (I4)
c     obj_a : Semi-major axis of object (n*R8)
c     obj_e : Eccentricity of object (n*R8)
c     obj_i : Inclination of object (n*R8)
c     obj_no: Longitude of node of object (n*R8)
c     obj_pe: Argument of perihelie of object (n*R8)
c     obj_m : Mean anomaly of object (n*R8)
c     obj_h : Absolute magnitude of object (n*R8)
c     obj_t : Time of elements of object (n*R8)
c     color : Array of colors (10*R8)
c     comp  : Dynamical componant the object belongs to (n*CH)
c     g     : Slope of object (R8)
c     ierr  : Error code (I4)
c                0 : nominal run, reached maximum number of objects
c               10 : unable to open filen
c               20 : error reading record
c               30 : end of file reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer*4
     $  n_max

      real*8
     $  obj_a(*), obj_e(*), obj_i(*), obj_no(*), obj_pe(*), obj_m(*),
     $  obj_h(*), obj_t(*), color(*),
     $  a, e, inc, node, peri, mt, h, g, jday, mag, r, alpha

      integer*4
     $  lun_d, n_obj, ierr

      character
     $  distri*(*), comp(*)*10, co*10

      save jday

      ierr = 0

c Hard coded slope for magnitude ! Bad boy !
      g = -0.12d0

c Open and read in object distribution
      n_obj = 0
 100  continue
         call read_obj (distri, lun_d, a, e, inc, mt,
     $     peri, node, h, jday, color, co, ierr)

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
         obj_a(n_obj) = a
         obj_e(n_obj) = e
         obj_i(n_obj) = inc
         obj_h(n_obj) = h
         obj_no(n_obj) = node
         obj_pe(n_obj) = peri
         obj_m(n_obj) = mt
         obj_t(n_obj) = jday
         if (co .ne. '          ') comp(n_obj) = co
         if (n_obj .ge. n_max) then
            return
         end if
         goto 100
 110  continue

      return
      end

      subroutine read_obj (filen, lun_in, a, e, i, capm,
     $  om, capom, h, jday, color, co, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine opens and reads in the object element file.
c Angles are returned in radian.
c Potentially use a common to return the H distribution parameters.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c Version 2 : For L7 data release, June 2010
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     filen : object element file name
c     lun_in: File unit
c
c OUTPUT
c     a     : Semi-major axis (R8)
c     e     : Eccentricity of orbit (R8)
c     i     : Inclination (R8)
c     capm  : Mean anomaly (R8)
c     om    : Argument of pericenter (R8)
c     capom : Longitude of node (R8)
c     h     : Absolute magnitude (R8)
c     jday  : Time of elements (R8)
c     color : Array of colors (10*R8)
c     co    : Dynamical componant the object belongs to (CH)
c     ierr  : Error code
c                0 : nominal run
c               10 : unable to open filen
c               20 : error reading record
c               30 : end of file reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      real*8
     $  a, e, i, capm, om, capom, h, Pi, drad, jday, color(*), jd

      integer*4
     $  nw_max

      parameter
     $  (Pi = 3.141592653589793238d0, drad = Pi/180.0D0, nw_max = 20)

      integer
     $  lun_in, ierr, j, nw, lw(nw_max)

      character
     $  line*100, filen*(*), word(nw_max)*80, co*10

      logical
     $  opened

      data opened /.false./

      save opened, jd

      ierr = 0
      if (.not. opened) then
         open (unit=lun_in, file=filen, status='old', err=1000)
         opened = .true.
         jd = -1.d0
         do j = 1, 10
            color(j) = 0.d0
         end do
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
      read (word(1), *) a
      read (word(2), *) e
      read (word(3), *) i
      i = i*drad
      read (word(4), *) capom
      read (word(5), *) om
      read (word(6), *) capm
      read (word(7), *) h
      capom = capom*drad
      om = om*drad
      capm = capm*drad
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

      end
