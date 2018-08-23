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
c 10 to 15. This range in reseved for them and won't be used by the
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

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c First define variables so they are accessible from a Python wrapper
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      block data modeldata

      real*8
     $  lambdaN, epoch_m,
     $  inamax, inamin, inimin, inimax

      common /com_inner/ inamax, inamin, inimin, inimax

      common /com_time/ epoch_m, lambdaN

      data
c One cannot use the form "sh /15.d0*drad/", so one needs to perform
c the computation before hand and give the result here, or do the
c multiplication in the program.
c
c Values for Inner
c
     $  inamax /39.d0/,
     $  inamin /37.d0/,
     $  inimin /7.d0/,
     $  inimax /20.d0/,
c
c Values of time and planet positions for all models
c
     $  lambdaN /5.489d0/,
     $  epoch_m /2453157.5d0/
      end

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c Now define working routines and functions
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
c This routine uses logical unit 10 to access the file containing the model.
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
     $  lun_m,                  ! Logical unit number, model description file
     $  lun_ll                  ! Logical unit number, logging
      real*8
     $  Pi,                     ! Pi
     $  TwoPi,                  ! 2*Pi
     $  drad                    ! Degree to radian convertion: Pi/180
c Set the values
      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0, lun_m = 20, lun_ll = 21)

c Internal storage
      integer*4
     $  flag,                   ! Tell invar_ecl_osc which direction to go
     $  i,                      ! Dummy index
     $  nparam                  ! Number of parameters for the function called
                                ! by routine incdism
      real*8
     $  beta_ah,                ! Slope of "a" distribution
     $  a0sl,                   ! Intermediate variable for "a" distribution
     $  a1sl,                   ! Intermediate variable for "a" distribution
     $  inamax,                 ! Lower "a" bound of inner belt
     $  inamin,                 ! Larger "a" bound of inner belt
     $  inimin,                 ! Lower limit of secular instability zone
     $  inimax,                 ! Upper limit of secular instability zone
     $  incmin,                 ! Lower limit of inclination distribution
     $  incmax,                 ! Upper limit of inclination distribution
     $  brown_params(5),        ! Parameters for the inclination distribution
     $  onecomp,                ! Returns probability of given inclination
     $  h_params(6),            ! Parameters for the H distribution
     $  size_dist_one,          ! Returns H according to exponential distrib
     $  q,                      ! Perihelion distance
     $  qmin,                   ! Lower limit of q distribution
     $  qmax,                   ! Upper limit of q distribution
     $  qhot,                   ! Returns the probability of given "q"
     $  param(10),              ! Temporary storage for distribution parameters
     $  epoch_m,                ! Epoch of elements [JD]
     $  lambdaN,                ! Longitude of Neptune at epoch
     $  random,                 ! Random number
     $  ran3,                   ! Random number generator
     $  r,                      ! Distance of object to Sun
     $  x,                      ! X cartesian coordinate of object
     $  y,                      ! Y cartesian coordinate of object
     $  z,                      ! Z cartesian coordinate of object
     $  color0(10)              ! Color parameters of model
      logical
     $  hot,                    ! Is this the hot component ?
     $  first,                  ! Tells if first call to routine
     $  log                     ! True if we want to log data

c Lightcurve and opposition surge effect parameters
      real*8
     $  gb0, ph0, period0, amp0

c Place some variables in common block so they can be accessed directly
c by a Python program.
      common /com_inner/ inamax, inamin, inimin, inimax
      common /com_time/ epoch_m, lambdaN

c Ran3, size_dist_one, onecomp and qhot are functions defined outside
c this routine
      external ran3, size_dist_one, onecomp, qhot

c Sets initial values
      data
     $  first /.true./,         ! First call
     $  gb0     /-0.12d0/,      ! Opposition surge effect
     $  ph0     / 0.00d0/,      ! Initial phase of lightcurve
     $  period0 / 0.60d0/,      ! Period of lightcurve
     $  amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

c Make sure variables are retained from call to call
      save

c This is the first call
      if (first) then
c Define the region of the inner belt
         incmin = 0.d0
         incmax = 180.d0*drad
         qmin = 34.d0
         qmax = 47.d0
c Define the "a" distribution.
         beta_ah = -2.5d0
         a0sl = inamin**(beta_ah + 1.d0)
         a1sl = inamax**(beta_ah + 1.d0)
c Reads in other parameters describing the model.
         open (unit=lun_m, file=filena, status='old', err=1000)
         read (lun_m, *) (h_params(i),i=1,3)
         read (lun_m, *) (h_params(i),i=4,6)
         read (lun_m, *) (brown_params(i),i=1,3)
         brown_params(2) = brown_params(2)*drad
         brown_params(3) = brown_params(3)*drad
         read (lun_m, *) (color0(i),i=1,10)
         read (lun_m, *) log
         close (lun_m)
c Writes a file describing the model that was used.
         open (unit=lun_ll, file='ModelUsed.dat', access='sequential',
     $     status='unknown')
         write (lun_ll, '(a)')
     $     'Inner belt model, version 1.0, 2013-05-14'
         write (lun_ll, '(''#'')')
         write (lun_ll, '(a,3(1x,f6.2))')
     $     '# H limits and slope for cold:', (h_params(i),i=1,3)
         write (lun_ll, '(a,3(1x,f6.2))')
     $     '# H limits and slope for hot: ', (h_params(i),i=4,6)
         write (lun_ll, '(a,3(1x,f6.2))')
     $     '# Width of cold, warm and hot:',
     $     (brown_params(i)/drad,i=1,3)
         write (lun_ll, '(''#'')')
         close (lun_ll)
c Change "first" so this is not called anymore
         first = .false.
      end if

 1100 continue
c
c Determination of "a" distribution, power-law
      random=ran3(seed)
      a = (a0sl + (a1sl-a0sl)*random)**(1.d0/(beta_ah+1.d0))
c
c Determination of "q" distribution
 1110 continue
      nparam = 0
      call incdism (seed, nparam, param, qmin, qmax, q,
     $  10, ierr, qhot)
c "q" can't exceed "a", can it ?
      if (q .gt. a) goto 1100
c Now get eccentricity
      e = 1.d0 - q/a
 1200 continue

c Select the component (hot or cold)
      random=ran3(seed)
      hot = (random .lt. brown_params(1))
c
c This is the hot component
      if (hot) then
         commen = 'h inner'
         nchar = 7
c
c Determination of "i" distribution
         param(1) = brown_params(3)
         nparam = 1
         call incdism (seed, nparam, param, incmin, incmax, inc,
     $     7, ierr, onecomp)
c
c This is the cold component
      else
         commen = 'c inner'
         nchar = 7
c
c Determination of "i" distribution
         param(1) = brown_params(2)
         nparam = 1
         call incdism (seed, nparam, param, incmin, incmax, inc,
     $     8, ierr, onecomp)
      end if

c Rejects if inclination in the secular instability zone
      if ((inc .gt. inimin*drad) .and. (inc .lt. inimax*drad))goto
     $  1200
c
c Angles: uniform distribution on allowable values
      random=ran3(seed)
      node = random*TwoPi
      random=ran3(seed)
      peri = random*TwoPi
      random=ran3(seed)
      m = random*TwoPi
c
c H-mag distribution: exponential law
      if (commen(1:1) .eq. 'h') then
         h = size_dist_one(seed, h_params(4))
      else
         h = size_dist_one(seed, h_params)
      end if

c Set up epoch for orbial elements
      epoch = epoch_m

c Define values for lightcurve and opposition surge effect
      gb = gb0
      ph = ph0
      period = period0
      amp = amp0

c Get colors for object
      do i = 1, 10
         color(i) = color0(i)
      end do

c The model above gives orbital elements with respect to the invariable
c plane reference frame (inclination 1$^\circ$ 35 13. 86 with respect to
c J2000 ecliptic plane with direction of ascending node at
c 107$^\circ$ 36 30. 8).
c The survey simulator expects the orbital elements with respect to the
c ecliptic, so convert them.
      flag = 1
      call invar_ecl_osc (flag, a, e, inc, node, peri, m,
     $  a, e, inc, node, peri, m, ierr)

c Store object if user requested
c Normally, we should explicitly open a file and write to its end, it
c seems like the pointer to the file is not retained from one call to
c the other, so simply use the default file assigned to the logical
c unit. In this case, the output file will be something like "fort.11"
      if (log) then
         call pos_cart(a,e,inc,node,peri,M,x,y,z)
         r = sqrt(x*x + y*y + z*z)
         open (unit=lun_ll, file='ModelUsed.dat', access='append',
     $     status='old')
         write(lun_ll,101) a,e,inc/drad,node/drad,peri/drad,M/drad,h,r,
     $     commen(1:nchar)
         close (lun_ll)
 101     format(6(f8.4,1x),f6.2,1x,f8.4,1x,a)
      end if

c Prepare return code
      ierr = 0

      return

 1000 continue
c If we get here, there is something really wrong, better return with
c panic code.
      ierr = -20
      return

      end

      real*8 function qhot (np, p, q)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c This function returns the probability of having perihelion distance "q"
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2013
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT
c     np    : Number of parameters describing the distribution function (I4)
c     p     : Array of parameters (np*R8)
c     q     : Perihelion distance (R8)
c
c OUTPUT
c     qhot  : Probability of q
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in) np
Cf2py intent(in) p
Cf2py intent(in) q
c
      implicit none

c Calling arguments
      integer np
      real*8 p(*), q

      qhot = 1./((1.+exp((35.-q)/0.5))*(1.+exp((q-40.)/0.5)))

      return
      end

      include 'ModelUtils.f'
