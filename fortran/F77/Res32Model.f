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
     $  lambdaN, epoch_m, emaxPLUT

      common /com_time/ epoch_m, lambdaN
      common /com_32/ emaxPLUT

      data
c One cannot use the form "sh /15.d0*drad/", so one needs to perform
c the computation before hand and give the result here, or do the
c multiplication in the program.
c
c Maximum eccentricity of plutinos
     $  emaxPLUT /0.32d0/,
c
c Values of time and planet positions for all models
c
     $  lambdaN /5.489d0/,
     $  epoch_m /2453157.5d0/
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
      real*8 a, e, inc, node, peri, M, epoch, h, color(10), gb, ph,
     $  period, amp
      character filena*(*), commen*100

c Some values better set up as parameters
      integer*4
     $  lun_m,                  ! Logical unit number for model description file
     $  lun_ll                  ! Logical unit number for logging
      real*8
     $  Pi,                     ! Pi
     $  TwoPi,                  ! 2*Pi
     $  drad                    ! Degree to radian convertion: Pi/180
c Set the values
      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0, lun_m = 20, lun_ll = 21)

c Internal storage
      character
     $  date*8,                 ! Date of execution
     $  time*10,                ! Time of execution
     $  zone*5                  ! Time zone
      integer*4
     $  values(8),              ! Values of date and time
     $  i,                      ! Dummy index
     $  j,                      ! Dummy index
     $  kozflag,                ! Flag for kozai resonance
     $  nparam                  ! Number of parameters for distribution
      real*8
     $  resamp,                 ! Resonance amplitude
     $  phi32,                  ! Resonant angle
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
     $  ecent,                  ! mean eccentricity
     $  ew,                     ! standard deviation of eccentricity
     $  emaxPLUT,               ! Maximum eccentricity of plutinos
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

c Lightcurve and opposition surge effect parameters
      real*8
     $  gb0, ph0, period0, amp0,
     $  color0(10)              ! Color parameters of model

c Vars for Kozai stuff 
      real*8
     $  kappa, etest, stest, hamtest
c, emin, emax, smin, smax
      real*8
     $  emina(78), emaxa(78), smina(78), smaxa(78), 
     $  wampa(78), hama(78), imina(78), imaxa(78)
      integer*4  nlines
c disturbing function constants (Wan and Huang 2007)
      real*8
     $  f1, f2, f3, f4, f7, f8, f20, f31, f53, f57,
     $  f33, f85, f34, f38, f87, f55, f94, f56, f97,
     $  f59, f60, f98

      parameter
     $  (f1=1.22929987494673, 
     $  f2=1.15279980000765, 
     $  f3=-4.61119920003061,
     $  f4=8.91299634372390, 
     $  f7=-101.627224295057,
     $  f8=99.3216246950413, 
     $  f20=34.9344527452960,
     $  f31=2.48400518330394,
     $  f53=8.26209324211937,
     $  f57=2.29309477909734,
     $  f33=1.30918774132387,
     $  f85=33.8114340769523,
     $  f34=-52.1519040380250,
     $  f38=3.33415568496546, 
     $  f87=22.2556656117726, 
     $  f55=-25.0365507216464,
     $  f94=151.734962414161, 
     $  f56=-302.645968671304,
     $  f97=3.84495442174114, 
     $  f59=166.317806052861, 
     $  f60=-94.5573022804512,
     $  f98=11.1996793328624)      

c Place some variables in common block so they can be accessed directly
c by a Python program.
      common /com_32/ emaxPLUT
      common /com_time/ epoch_m, lambdaN

c Ran3, size_dist_one and gasdev are functions defined outside
c this routine.
      external ran3, size_dist_one, gasdev, onecomp

c Sets initial values
      data 
     $  first   /.true./,         ! First call
     $  gb0     /-0.12d0/,      ! Opposition surge effect
     $  ph0     / 0.00d0/,      ! Initial phase of lightcurve
     $  period0 / 0.60d0/,      ! Period of lightcurve
     $  amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

      save

c This is the first call
      if (first) then
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
         read(lun_m,*) ew
         read(lun_m,*) ecent
         read(lun_m,*) libampc
         read(lun_m,*) fkoz
         read(lun_m,*) nlines
         read(lun_m,*) kappa
         do j = 1, nlines
            read(lun_m,*) emina(j), emaxa(j), imina(j), imaxa(j),
     $        hama(j), wampa(j)
            smina(j)=sin(imina(j)/(2.0d0*drad))
            smaxa(j)=sin(imaxa(j)/(2.0d0*drad))
            wampa(j)=wampa(j)*drad
         end do
         read (lun_m, *) (color0(i),i=1,10)
         read (lun_m, *) log
         close(lun_m)

c Writes a file describing the model that was used.
         open (unit=lun_ll, file='ModelUsed', access='sequential',
     $     status='unknown')
         write (lun_ll, '(a)')
     $     '# 3:2 resonance model, version 1.0, 2013-06-12'
         call date_and_time(date, time, zone, values)
         write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ',
     $     date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' 
     $     //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
         write (lun_ll, '(''#'')')
         write (lun_ll, '(a,5(1x,f6.2))')
     $        '# H limits and slopes:', (h_params(i),i=1,5)
         write (lun_ll, '(a,1(1x,f6.2),a)')
     $        '# Width of inc:              ', sg2deg
         write (lun_ll, '(''#'')')
         close (lun_ll)
c Change "first" so this is not called anymore
         first = .false.
      end if

c decide if it's a Kozai Plutino
      random=ran3(seed)
      if (random .lt. fkoz) goto 3241
c jump to Kozai section, otherwise continue with normal Plutino

      kozflag=0
c NON-KOZAI PLUTINO SECTION
 1051 continue
c ADDED GASDEV
      random=ran3(seed)
c nominal is 0.18, 0.06
      e = gasdev(ecent,ew,seed)
      if (e .lt. 0.0d0) goto 1051

c pick 'a' and re-draw if a/e outside of bounds
      random=ran3(seed)
c center and shape from stability plots in Tiscareno paper
      a= 39.45d0 + (random-0.5d0)*0.4d0
      if (a .gt. (39.45d0 + 4.0d0/3.0d0*(e - 0.01d0)) ) goto 1051
      if (a .lt. (39.45d0 - 4.0d0/3.0d0*(e - 0.01d0)) ) goto 1051

c        cutoff if Uranus-approaching.  
      q = a*(1.0d0-e)
      if (q .lt. 22.d0) goto 1051

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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             Need random resonant argument, with libration amplitude
c             phi32 = 3*lambda - 2*lambdaN - longperi
c             phi32 = 3*(M+peri+node) - 2*lambdaN - (node+peri)
c             Here use Pluto's libration amplitude  (85 degrees)
cccc          phi32 = Pi + 2.0d0*(random - 0.5d0)*85.d0*drad

c try a triangle (JMP) based on Lykawka's histogram.
c
 1054 continue
      resamp = ran3(seed)*110.0d0+20.0d0
      random = ran3(seed)
      if (resamp .lt. libampc) then
         if (random .gt. resamp/libampc) goto 1054
      else
         if (random .gt. (130.d0-resamp)/(130.d0-libampc)) goto 1054
      end if

ccccccccccccccccccccccc Now choose phi32

      random=ran3(seed)
      phi32 = Pi + sin(2.0d0*Pi*random)*resamp*drad

      random=ran3(seed)
      if (random .gt. 0.5) phi32 = phi32 - TwoPi ! Object is trailing

c             Set argument of pericenter based on resonant angle
      peri = 0.5d0*(phi32 - 3.d0*m +2.d0*lambdaN)  - node
      call ztopi(peri)

      goto 1999

c ----------------------------------------------
c KOZAI
c ----------------------------------------------
 3241 continue
      kozflag=1

c randomly pick a line from the koztable file
      random=ran3(seed)*float(nlines)
      j=aint(random)

c pick peri from 90-omegamax, effectively picking Hamiltonian
      random=ran3(seed)
      peri=Pi/2.0d0 + (Pi/2.0d0 - wampa(j))*sin(2.0d0*Pi*random)

c pick e using emin, emax, omega (peri) and hamiltonian
c start half from "outside" hamiltonian curve, half "inside"
      random=ran3(seed)
      if (random .lt. 0.5) then 
         etest=0.32d0
 4486    stest=sqrt(0.5d0*(1.0d0 - 
     $     sqrt(kappa**2.0d0 / (1.0d0 - etest**2.0d0))))
         hamtest = f1 - 
     $     f31*etest + 
     $     (f2+f53)*etest**2.0d0 + 
     $     f3*stest**2.0d0 - 
     $     (f33+f85)*etest**3.0d0 - 
     $     f34*etest*stest**2.0d0 + 
     $     (f4+f94+f55)*etest**4.0d0 + 
     $     (f7+f56)*etest**2.0*stest**2.0d0 + 
     $     f8*stest**4.0d0 + 
     $     (f57*stest**2.0d0 - 
     $     (f38+f87)*etest*stest**2.0d0 + 
     $     (f20+f97+f59)*etest**2.0d0*stest**2.0d0 + 
     $     f60*stest**4.0d0)*cos(2.0d0*peri) + 
     $     f98*stest**4.0d0*cos(4.0d0*peri)
c        write(6,8787) etest,emaxa(j),hamtest,hama(j),peri,wampa(j)
c 8787  format (6(f20.17,1x))
         if ((hamtest .gt. hama(j)) .and. (etest .gt. 0.15d0)) then
            etest=etest-0.001d0
            if (etest .le. 0.15d0) goto 3241
            goto 4486
         end if 
      else
         etest=0.15d0
 4487    stest=sqrt(0.5d0*(1.0d0 - 
     $     sqrt(kappa**2.0d0 / (1.0d0 - etest**2.0d0))))
         hamtest = f1 - 
     $     f31*etest + 
     $     (f2+f53)*etest**2.0d0 + 
     $     f3*stest**2.0d0 - 
     $     (f33+f85)*etest**3.0d0 - 
     $     f34*etest*stest**2.0d0 + 
     $     (f4+f94+f55)*etest**4.0d0 + 
     $     (f7+f56)*etest**2.0*stest**2.0d0 + 
     $     f8*stest**4.0d0 + 
     $     (f57*stest**2.0d0 - 
     $     (f38+f87)*etest*stest**2.0d0 + 
     $     (f20+f97+f59)*etest**2.0d0*stest**2.0d0 + 
     $     f60*stest**4.0d0)*cos(2.0d0*peri) + 
     $     f98*stest**4.0d0*cos(4.0d0*peri)
c        write(6,8788) etest,emaxa(j),hamtest,hama(j),peri,wampa(j)
c 8788   format (6(f20.17,1x))
         if ((hamtest .gt. hama(j)) .and. (etest .lt. 0.32d0)) then
            etest=etest+0.001d0
            if (etest .ge. 0.32d0) goto 3241
            goto 4487
         end if 
      end if
      e=etest

c inc is in radians here
      inc= 2.0d0*asin(stest)

c flip half to 360 - peri 
      random=ran3(seed)
      if (random .lt. 0.5) peri=2.0d0*Pi - peri 

c pick 'a' and re-draw if a/e outside of bounds
      random=ran3(seed)
c center and shape from stability plots in Tiscareno paper
      a= 39.45 + (random-0.5d0)*0.4d0
      if (a .gt. (39.45d0 + 4.0d0/3.0d0*(e - 0.01d0)) ) goto 3241
      if (a .lt. (39.45d0 - 4.0d0/3.0d0*(e - 0.01d0)) ) goto 3241

c M picked randomly.  
      random=ran3(seed)
      m = random*TwoPi

c Libration, should be same as for non-Kozai Plutino
c try a triangle (JMP) based on Lykawka's histogram.
 1055 continue
      resamp = ran3(seed)*110.0d0+20.0d0
      random = ran3(seed)
      if (resamp .lt. libampc) then
         if (random .gt. resamp/libampc) goto 1055
      else
         if (random .gt. (130.d0-resamp)/(130.d0-libampc)) goto 1055
      end if

      random=ran3(seed)
      phi32 = Pi + sin(2.0d0*Pi*random)*resamp*drad

      random=ran3(seed)
      if (random .gt. 0.5) phi32 = phi32 - TwoPi ! Object is trailing

c             Set node based on resonant angle
      node = 0.5d0*(phi32 - 3.d0*m) + lambdaN  - peri
      call ztopi(node)
      call ztopi(peri)

      goto 1999

 1999 continue

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

c Define values for lightcurve and opposition surge effect
      gb = gb0
      ph = ph0
      period = period0
      amp = amp0

c Get colors for object
      do i = 1, 10
         color(i) = color0(i)
      end do

c Stores informations on resonance
      commen = '03:02 '
      write (commen(7:17), '(f10.6,1x)') phi32
      write (commen(18:28), '(f10.5,1x)') resamp
      write (commen(29:29), '(i1)') kozflag
      nchar = 29

c Stores object if user requested
      if (log) then
         call pos_cart(a,e,inc,node,peri,M,x,y,z)
         r = sqrt(x*x + y*y + z*z)
         open (unit=lun_ll, file='ModelUsed.dat', access='append',
     $     status='old')
         write(lun_ll,101) a,e,inc/drad,node/drad,peri/drad,M/drad,h,r,
     $     x, y, z, commen(1:nchar)
         close (lun_ll)
 101     format(6(f8.4,1x),f6.2,1x,f8.4,3(1x,f7.3),1x,a)
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
