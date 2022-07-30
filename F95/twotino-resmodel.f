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
     $  lambdaN, epoch_m, emax52

      common /com_time/ epoch_m, lambdaN
      common /com_52/ emax52

      data
c One cannot use the form "sh /15.d0*drad/", so one needs to perform
c the computation before hand and give the result here, or do the
c multiplication in the program.
c
c Maximum eccentricity of plutinos
     $  emax52 /0.32d0/,
c
c Values of time and planet positions for all models
c
     $  lambdaN /5.489d0/,
c for epoch_m /2456839.5d0/
     $  lambdaN /5.876d0/,
c     $  lambdaN /5.83917354547d0/,
c     $  epoch_m /2456505.5d0/
     $  epoch_m /2453157.5d0/
c 2014-Jul-01 2456839.5, lambdaN /5.876d0/,
     $  epoch_m /2456839.5d0/
      end

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c Now define working routines and functions
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

c  *********************************************************************
c  *********************************************************************
c  Subroutine res_3_2
c
c  currently generates objects in the 3:2 mean motion resonance
c  Based on what was in plut.f from kozai_plutinos
c
c  Set up to generate a from a uniform +/- about the res center
c  e from a gaussian distribution (truncated at e=0)
c  i from a sin(i)*gaussian distribution
c  resonant amplitude from a triangle shaped thing
c
c  Should be replaced/modified as the user wishes
c      
c  Kat Volk May 2013
c
c  2013-05-14: J.-M. Petit adapted to new GiMeObj routine
c
c  2013-07-01: J.-M. Petit modified API to fit new Drivers.f
c
c  2013-09-20: J.-M. Petit changed H distributino to include divot
c
c  *********************************************************************
c  *********************************************************************

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
      real*8 a, e, inc, node, peri, M, epoch, h, color(*), gb, ph,
     $  period, amp
      character filena*(*), commen*100, log_file*20
      parameter (log_file = 'ModelUsed')
c Some values better set up as parameters
      integer*4
     $  lun_m,                  ! Logical unit number for model description file
     $  lun_l,                  ! Logical unit number for logging
     $  kozflag
      real*8
     $  Pi,                     ! Pi
     $  TwoPi,                  ! 2*Pi
     $  drad                    ! Degree to radian convertion: Pi/180
c Set the values
      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi, drad = Pi/180.0d0,
     $  lun_m = 90, lun_l = 95)

c Internal storage
      integer*4
     $  i,                      ! Dummy index
     $  j,                       ! Dummy index
     $ nparam,
     $ lts                      !sym = 0, leading = 1, trailing = 2
      real*8
     $  resamp,                 ! Resonance amplitude
     $  fsym,                   ! symmetric fraction
     $  phi52,                  ! Resonant angle
     $  size_dist_one,          ! Returns H according to exponential distrib
     $  h_params(6),            ! Parameters for the H distribution
     $  params(6),              ! Parameters for some distribution
     $  incmin,                 ! Lower limit of inclination distribution
     $  incmax,                 ! Upper limit of inclination distribution
     $  epoch_m,                ! Epoch of elements [JD]
     $  random,                 ! Random number
     $  ran3,                   ! Random number generator
     $  gasdev,                 ! Returns a random number with gaussian
                                ! distribution of given width and center
     $  q,                      ! Perihelion distance
     $  lambdaN,                ! Neptune's mean longitude at some reference epoch
     $  emin,                  ! mean eccentricity
     $  emax,                     ! standard deviation of eccentricity
     $  emins,                  ! mean eccentricity
     $  emaxs,                     ! standard deviation of eccentricity
     $  emax52,               ! Maximum eccentricity of plutinos
     $  fleading,               ! fraction of asymm lib that are in leading island
     $  sg2deg,                 ! standard deviation in inclination (degree)
     $  sg2,                    ! standard deviation in inclination (radian)
     $  libampc,                ! center of the libration amp dist
     $  fkoz,                   ! Fraction of kozai resonators
     $  h0,                     ! lowest value for h distribution
     $  h1,                     ! knee value for h distribution
     $  h2,                     ! highest value for h distribution
     $  slope,                  ! slope of size distribution, before divot
     $  slope2,                 ! slope of size distribution, past divot
     $  ct,                     ! constrast factor at divot
     $  ckc,                    ! nb/(nb+ns)
     $  nb,                     ! number of ojbect brighter than divot
     $  ns,                     ! number of ojbect fainter than divot
     $  h0s10,                  ! normalizing factor at low end of H-distrib
     $  h1s10,                  ! normalizing factor at high end of H-distrib
     $  r,                      ! Distance of object to Sun
     $  x,                      ! X cartesian coordinate of object
     $  y,                      ! Y cartesian coordinate of object
     $  z,                      ! Z cartesian coordinate of object
     $  onecomp                 ! Name of function for inclination distribution
      logical
     $  first,                  ! Tells if first call to routine
     $  log                     ! True if we want to log data  

      real*8 libc, ymax, ampmax, mslope, binter, temp, temp2, libcl,
     $   libcf, amplowest, ldl
c Place some variables in common block so they can be accessed directly
c by a Python program.
      common /com_52/ emax52
      common /com_time/ epoch_m, lambdaN

c Ran3, size_dist_one and gasdev are functions defined outside
c this routine.
      external ran3, size_dist_one, gasdev, onecomp

c Sets initial values
      data first /.true./

      save

c This is the first call
      if (first) then
c Define values for lightcurve and opposition surge effect
         gb = -0.12d0
         ph = 0.d0
         period = 0.6d0
         amp = 0.d0
c Reads in other parameters describing the model.
         open(unit=lun_m, file=filena, status="unknown", err=1000)
         read(lun_m,*) sg2deg
         sg2=sg2deg*drad
         read(lun_m,*) h0
         read(lun_m,*) h1
         read(lun_m,*) h2
         read(lun_m,*) slope
         read(lun_m,*) slope2
         h_params(1) = h0
         h_params(2) = h1
         h_params(3) = h2
         h_params(4) = slope
         h_params(5) = slope2
         read(lun_m,*) ct
         nb = (10.d0**(slope*h1) - 10.d0**(slope*h0))/slope
         ns =10.d0**((slope-slope2)*h1)/slope2
     $     *(10.d0**(slope2*h2) - 10.d0**(slope2*h1))/ct
         ckc = nb/(nb + ns)
         read(lun_m,*) emin
         read(lun_m,*) emax
         read(lun_m,*) emins
         read(lun_m,*) emaxs
         read(lun_m,*) fleading
         read(lun_m,*) fsym
         read (lun_m, *) (color(i),i=1,10)
         read (lun_m, *) log
         close(lun_m)

c Writes a file describing the model that was used.
c         open (unit=lun_m, file='ModelUsed', status='unknown')
c         write (lun_m, '(a)')
c     $     '3:2 resonance model, version 1.0, 2013-06-12'
c         close (lun_m)
c Change "first" so this is not called anymore

         first = .false.
      end if


c pick 'a' and re-draw if a/e outside of bounds
         random=ran3(seed)
         a= 47.8 + (random-0.5)*2.0 * 0.20


c Pull an inclination
      random = ran3(seed)
c J-M's code:
      params(1) = sg2
      nparam = 1
      incmin = (0.0d0*drad)
      incmax = (90.0d0*drad)
      call incdism (seed, nparam, params, incmin, incmax, inc,
     $  1, ierr, onecomp)

c node and M picked randomly.  
      random=ran3(seed)
      node = random*TwoPi
      random=ran3(seed)
      m = random*TwoPi


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c have to pick whether it's symmetric or asymmetric, and which island
c it orbits

         random=ran3(seed)
         if (random .le. fsym) goto 1235 
c symmetric, otherwise asymmetric and code continues
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


1051     continue

c new asymmetric rules:
c pick an e.
c maximum libration amplitude for that given e is roughly:
c 81.6 - 1.5631/sinh(e/1.636) based solely off a rough fit to the figure
c in Nesvorny & Roig 2001
c amp;itude is then uniformly selected between 0 and that value
c based off the same figure, the libration center is given by:
c 72.75 - 33.1*x - 15.32*x**2 + 3.1616/x

c pick e
         e = ran3(seed)*(emax-emin)+emin
c give a slope distribution
133      e =  gasdev(0.275d0,0.060d0,seed)
         if (e .lt. emin) goto 133
c        cutoff if Uranus-approaching. , reselect e 

         q = a*(1.0d0-e)
         if (q .lt. 22.000) goto 1051         

c pick a libration center


         
520         if(e .lt. 0.05) then
            libc = 140d0 + (1d0-ran3(seed))*10d0
            ampmax = 15d0
         else
c            libc = 72.75d0 - 33.1d0*e - 15.32*e*e + 3.1616/e
c            libc = 57.391d0 +69.9391d0*e -158.766d0*e*e + 3.85541d0/e
c        get the value from the line of e v.s. libc
c            libcl = 174.986d0 -551.752d0*e+792.648d0*e*e -2.77753d0/e 
            libcl = 123.876d0 -272.425d0*e+325.087d0*e*e +0.0799027d0/e
            libcl = 118.309d0 -243.407d0*e+279.281d0*e*e +0.40247d0/e
            libcl = 133.157d0 -315.272d0*e+377.082d0*e*e -0.491667/e
c        get the libc
            random=ran3(seed)
            libc = libcl+ (1-sin(0.5d0*Pi*random))*20

c            ampmax = 81.6d0 - 1.5631d0/sinh(e/1.636d0)
c            ampmax = -8.85454d0 -1.19166d0/sinh(e/(-5.48914d0))
c            ampmax = 689.264d0 -3.40183d0*libc+0.00302945d0*libc*libc
c     $         -31567.5d0/libc
            ampmax = -403.632d0 +9.09917d0*libc-0.0442498d0*libc*libc
     $         -0.0883975d0/libc
            amplowest = 2006.19d0 -27.505d0*libc+0.125945d0*libc*libc
     $         -48563.1d0/libc
            amplowest = 2598.64d0 -176.863d0*libc+4.7044d0*libc**2 
     $        -0.0611052d0*libc**3 +0.000387165d0*libc**4
     $        -9.53767e-07*libc**5
            amplowest = 79.031d0*exp(-(libc-121.3435d0)**2/
     $         (2*15.51349d0**2))
            lts = 1
            random = ran3(seed)
            if (random .gt. fleading) then
                 libc=360.0d0-libc
                 lts = 2 !trailling
            endif

         end if

1054     continue

c zzzz
c         resamp = ran3(seed)*ampmax 
c         resamp = ampmax + (sin(2.0d0*Pi*ran3(seed))*20 -10)

C	pick which side of the distribution center to be on
c         temp = ran3(seed)
c         if(temp .gt. 0.5) resamp = ampmax+5d0-resamp
c         resamp = ran3(seed)*20 -10 + ampmax
c         resamp = resamp+ ran3(seed)*20 -10
         resamp = ampmax + (ran3(seed)-1)*(ampmax - amplowest) 
         resamp = ampmax - sin(0.5d0*Pi*ran3(seed))*ran3(seed)*
     $      (ampmax - amplowest)
222      random =  ran3(seed)
         ldl = 1.0d0 - 0.75d0*random
         if (ran3(seed) .gt. ldl) goto 222
         resamp = ampmax - (random)*(ampmax - amplowest)

            
c         resamp = ampmax - (ampmax - amplowest)*
c     $          (1-sin(0.5d0*Pi*ran3(seed)))
c     $          (cos(Pi/3.0d0 +0.16666d0*Pi*ran3(seed))+0.5)
c         if (resamp .lt. 0.000) goto 1051 
c pick whether leading or trailing
c         lts = 1
c         random = ran3(seed)
c          if (random .gt. fleading) then
c                 libc=360.0d0-libc
c                 lts = 2 !trailling
c          endif


ccccccccccccccccccccccc Now choose phi

c         random=ran3(seed)
c         phi52 = libc*drad + sin(2.0d0*Pi*random)*resamp*drad 
c         libc2 = libc+ (1-sin(0.5d0*Pi*random))*40-20
c         if(libc .gt. 180.0d0) libc2 = libc+ (sin(0.5d0*Pi
c     $     *random))*30-15

         phi52 = (libc + 2*resamp*(ran3(seed) - 0.5))*drad
c     $     *15*drad*ran3(seed)
c             Set argument of pericenter based on resonant angle
         peri = phi52 - 2.0d0*m + lambdaN  - node
         libcf = libc
         call zero2pi(peri)
         goto 1999

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c symmetric
1235     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c symmetric, libration centers are all 180
c            pick libration amplitudes uniformly from 125-165
c            pick e uniformly from 0.05 -> 0.3 (Chiang&Jordan02)
         lts = 0

c pick a libration center
         libampc=180.0d0

1052     continue
c pick e
         e = ran3(seed)*(emaxs-emins)+emins
c        cutoff if Uranus-approaching.  
         q = a*(1.0d0-e)
         if (q .lt. 22.000) goto 1052       

c pick a libration amplitude
          resamp = ran3(seed)*27d0 + 140d0

ccccccccccccccccccccccc Now choose phi

         random=ran3(seed)
c         phi52 = libampc*drad + sin(2.0d0*Pi*random)*resamp*drad
         phi52 = libampc*drad + (2*random-1)*resamp*drad
c             Set argument of pericenter based on resonant angle
         peri = phi52 - 2.0d0*m + lambdaN  - node
         libcf = libampc
         call zero2pi(peri)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c -----------------------------------------------
 1999    continue


c H-mag distribution
      random = ran3(seed)
      if (random .lt. ckc) then
         slope = h_params(4)
         h0 = h_params(1)
         h1 = h_params(2)
      else
         slope = h_params(5)
         h0 = h_params(2)
         h1 = h_params(3)
      end if
      h0s10 = 10.d0**(slope*h0)
      h1s10 = 10.d0**(slope*h1)
      random = ran3(seed)
      h = log10(random*(h1s10 - h0s10) + h0s10)/slope

      epoch = epoch_m

c Stores informations on resonance
      commen = '02:01 '
      kozflag = 0
      write (commen(7:17), '(f10.6,1x)') phi52
      write (commen(18:28), '(f10.5,1x)') resamp
      write (commen(29:29), '(i1)') lts
      write (commen(31:41), '(f10.6,1x)') libcf
      nchar = 41

c Stores object if user requested
      if (log) then
         call pos_cart(a,e,inc,node,peri,M,x,y,z)
         r = sqrt(x*x + y*y + z*z)
         open(unit=lun_l, file=log_file, position='append')
         write(lun_l,101) a,e,inc/drad,node/drad,peri/drad,M/drad,h,r,
     $     x, y, z, commen(1:nchar)
 101     format(6(f8.4,1x),f6.2,1x,f8.4,3(1x,f7.3),1x,a)
         close(unit=lun_l)
      end if
      ierr = 0

      return

 1000 continue
      ierr = -20
      return

      end

c--------------------------------------------
c   gasdev
c********************************************************************
c  Given a center, width and seed return a value drawn from a gaussian

      FUNCTION gasdev(x0,sigma,rs)

      implicit none
      real*8 x0, sigma
      integer*4 iset, rs
      real*8 v1, v2, rsq, gasdev, gset
      real*8 fac, random, ran3
      external ran3
      SAVE iset,gset
      
      if  (iset.eq.0) then
 12      continue
             random = ran3(rs)
             v1= 2d0*random - 1d0
             random = ran3(rs)
             v2= 2d0*random - 1d0
             rsq=v1*v1+v2*v2
         if (rsq.ge.1.0) goto 12 

         fac=sqrt(-2d0*log(rsq)/rsq)
         gset=v1*fac
         iset=1
         gasdev=v2*fac*sigma+x0
         return 
      else
         iset=0
         gasdev=gset*sigma+x0
         return
      endif
      END

      include 'ModelUtils.f'
















