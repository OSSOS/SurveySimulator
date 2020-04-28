module gimeobjut

  use datadec
  use elemutils
  use rot
  use modelutils
  use ioutils

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! First define variables so they are accessible from a Python wrapper
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  real (kind=8) :: lambdaN, epoch_m, inamax, inamin, inimin, inimax
  common /com_inner/ inamax, inamin, inimin, inimax
  common /com_time/ epoch_m, lambdaN

! Values for Inner
!
  data inamax /39.d0/, inamin /37.d0/, inimin /7.d0/, inimax /20.d0/, &
!
! Values of time and planet positions for all models
!
       lambdaN /5.489d0/, epoch_m /2453157.5d0/

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
! 10 to 15. This range in reseved for them and won't be used by the
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
! File generated on 2020-04-16
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
! Version 2 : April 2020
!             Upgrade to Fortran 95 syntax
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
! This routine uses logical unit 10 to access the file containing the model.
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
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(out) :: ierr, nchar
    type(t_orb_m), intent(out) :: o_m
    real (kind=8), intent(out) :: epoch, h, color(10), gb, ph, period, amp
    character(*), intent(in) :: filena
    character(100), intent(out) :: commen

! Some values better set up as parameters
    type(t_v3d) :: p
    type(func_holder) :: fp
    integer, parameter :: &
         lun_m = 20,            &! Logical unit number for data file reading
         lun_ll = 21             ! Logical unit number for logging
    real (kind=8), parameter :: &
         Pi = 3.141592653589793238d0, &! Pi
         TwoPi = 2.0d0*Pi,      &! 2*Pi
         drad = Pi/180.0d0       ! Degree to radian convertion: Pi/180

! Internal storage
    character(5) :: zone         ! Time zone
    character(8) :: date         ! Date of execution
    character(10) :: time        ! Time of execution
    integer (kind=4), save :: &
         values(8),             &! Date and time of execution
         flag,                  &! Tell invar_ecl_osc which direction to go
         i,                     &! Dummy index
         nparam                  ! Number of parameters for the function called
                                 ! by routine incdism
    real (kind=8), save :: &
         brown_params(5),       &! Parameters for the inclination distribution
         h_params(6),           &! Parameters for the H distribution
         q,                     &! Perihelion distance
         qmin,                  &! Lower limit of q distribution
         qmax,                  &! Upper limit of q distribution
         param(10),             &! Temporary storage for distribution parameters
         random,                &! Random number
         r,                     &! Distance of object to Sun
         color0(10)              ! Color parameters of model
    logical, save :: &
         hot,                   &! Is this the hot component ?
         first,                 &! Tells if first call to routine
         log                     ! True if we want to log data
    real (kind=8) :: &
         inamax,                &! Lower "a" bound of inner belt
         inamin,                &! Larger "a" bound of inner belt
         inimin,                &! Lower limit of secular instability zone
         inimax,                &! Upper limit of secular instability zone
         incmin,                &! Lower limit of inclination distribution
         incmax,                &! Upper limit of inclination distribution
         epoch_m,               &! Epoch of elements [JD]
         lambdaN                 ! Longitude of Neptune at epoch
! Lightcurve and opposition surge effect parameters
    real (kind=8), save :: gb0, ph0, period0, amp0

! Place some variables in common block so they can be accessed directly
! by a Python program.
    common /com_inner/ inamax, inamin, inimin, inimax
    common /com_time/ epoch_m, lambdaN

! Sets initial values
    data &
         first /.true./,        &! First call
         gb0     /-0.12d0/,     &! Opposition surge effect
         ph0     / 0.00d0/,     &! Initial phase of lightcurve
         period0 / 0.60d0/,     &! Period of lightcurve
         amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

! This is the first call
    if (first) then
! Define the region of the inner belt
       inamin = 37.d0         ! lowest  value of a (AU)
       inamax = 39.d0         ! largest value of a
       qmin   = 34.d0
       qmax   = 39.d0
       incmin = 0.d0  *drad
       incmax = 180.d0*drad
       inimin = 7.d0  *drad   ! This cuts out 7-20 deg inclinations (nu 8)
       inimax = 20.d0 *drad
! Reads in other parameters describing the model.
       open (unit=lun_m, file=filena, status='old', err=1000)
       read (lun_m, *) (h_params(i),i=1,3)     ! Size distribution parameters
       read (lun_m, *) (brown_params(i),i=1,1) ! Read i distribution width 
       brown_params(1) = brown_params(1)*drad  ! convert to radians
       read (lun_m, *) (color0(i),i=1,10)      ! read color array
       read (lun_m, *) log          ! logical variable, turn on drawing log
       close (lun_m)
! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a)') &
            'Inner belt model, version 1.0, 2013-05-14'
       call date_and_time(date, time, zone, values)
       write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ', &
            date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' &
            //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# H limits and slope for cold:', (h_params(i),i=1,3)
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# H limits and slope for hot: ', (h_params(i),i=4,6)
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# Width of cold, warm and hot:', &
            (brown_params(i)/drad,i=1,3)
       write (lun_ll, '(''#'')')
       close (lun_ll)
! Change "first" so this is not called anymore
       first = .false.
    end if

1100 continue
!
! Determination of "a" distribution, uniform from a_min to a_max
    random=ran_3(seed)
    o_m%a = inamin + (inamax-inamin)*random
!
! Determination of simple "q" distribution. This is NOT a good algorithm.
1110 continue
    random=ran_3(seed)
    q = qmin + (qmax-qmin)*random
    if (q .gt. o_m%a) goto 1100 ! q can't exceed a. Redraw another a if so
    o_m%e = 1.d0 - q/o_m%a

1200 continue
! Determination of "i" distribution. Width was read as brown_params(1)
    param(1) = brown_params(1)
    nparam = 1
! This utitily returns, when there is one parameter, an inc drawn from
! a sin(i)*gaussian inclination distribution.
! The 7 below is a unique code assigned to this distribution for speed.
    fp%f_ptr => onecomp
    call incdism (seed, nparam, param, incmin, incmax, o_m%inc, &
         7, ierr, fp)
    commen = 'h inner'
    nchar = 7
!
! Rejects if inclination in the secular instability zone
    if ((o_m%inc .gt. inimin*drad) .and. (o_m%inc .lt. inimax*drad)) goto 1200
!
! Angles: uniform distribution on allowable values
    random=ran_3(seed)
    o_m%node = random*TwoPi
    random=ran_3(seed)
    o_m%peri = random*TwoPi
    random=ran_3(seed)
    o_m%m = random*TwoPi
!
! H-mag distribution: exponential law
    h = size_dist_one(seed, h_params)

! Set up epoch for orbial elements
    epoch = epoch_m

! Define values for lightcurve and opposition surge effect
    gb = gb0
    ph = ph0
    period = period0
    amp = amp0

! Get colors for object
    do i = 1, 10
       color(i) = color0(i)
    end do

! The model above gives orbital elements with respect to the invariable
! plane reference frame (inclination 1$^\circ$ 35 13. 86 with respect to
! J2000 ecliptic plane with direction of ascending node at
! 107$^\circ$ 36 30. 8).
! The survey simulator expects the orbital elements with respect to the
! ecliptic, so convert them.
    flag = 1
    call invar_ecl_osc (flag, o_m, o_m, ierr)

! Store object if user requested
! Normally, we should explicitly open a file and write to its end, it
! seems like the pointer to the file is not retained from one call to
! the other, so simply use the default file assigned to the logical
! unit. In this case, the output file will be something like "fort.11"
    if (log) then
       call pos_cart(o_m, p)
       r = sqrt(p%x*p%x + p%y*p%y + p%z*p%z)
       open (unit=lun_ll, file='ModelUsed.dat', access='append', &
            status='old')
       write(lun_ll,101) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
            o_m%peri/drad, o_m%m/drad, h, r, commen(1:nchar)
       close (lun_ll)
101    format(6(f8.4,1x),f6.2,1x,f8.4,1x,a)
    end if

! Prepare return code
    ierr = 0

    return

1000 continue
! If we get here, there is something really wrong, better return with
! panic code.
    ierr = -20
    return

  end subroutine GiMeObj

end module gimeobjut
