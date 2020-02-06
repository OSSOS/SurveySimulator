c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Survey Simulator set of routines, version 2.0 for OSSOS
c
c Based on SurveySim.f from CFEPS, as of 2013-01-15
c
c Improvements include use of new efficiency and pointings file format,
c various efficiency functions for a given pointing as a function of rate of
c motion and colors for simulated objects such that surveys have not been
c shifted in magnitude to a common filter.
c
c 2013-04-26: Modified GetSurvey, read_sur and read_eff to read in new
c             efficiency and pointings files.
c             Modified eta_raw to support multiple efficiency functions as
c             functions of rate of motion.
c             Modified detect1 to add a 'color' array in the argument list to
c             account for colors of objects and various filters used in the
c             surveys, and change the sequence of rejections to:
c               - Check app. mag. against limiting mag.
c               - Check FOV
c               - Check for filling factor
c               - Check rate of motion
c               - Check efficiency on randomized mag.
c               - Check for tracking
c
c 2013-05-03: Added two more arguments at end of list of detect1: 'ic' for the
c             index of the color of the detection filter, and 'surna', a
c             3-character string with the name of the survey that detected the
c             object.
c             Updated content of efficiency files for CFEPS to eliminate
c             remaining inconsistancies.
c
c 2013-05-08: Changed name of entry subroutine from detect1 to detos1.
c
c 2013-07-01: Modified API for detos1 to return a single COMMENTS string
c             argument and an integer NCHAR specifying how many characeter to
c             output.
c
c 2013-09-18: Changed API of detos1 to return intrinsic mag and averaged
c             randomized mag instead of randomized mag and dmag.
c             Changed algorithm to check efficiency based on intrinsic
c             mag, using full efficiency function down to 1%, and finally
c             reject object if averaged randomized mag is fainter than
c             40% limit.
c             Discarded use of eta_trust.
c
c 2015-01-12: Changed API of detos1 to return surmized absolute magnitude
c             obtained from averaged randomized mag. Hx and m_int are returned
c             in the band filter used for the H distribution (X), while m_rand
c             and H_rand are returned in the band filter of the survey the
c             object is detected in.
c
c 2016-03-16: Changed detos1 to use 'poly', 'ears' or 'rect' labels instead of
c             the old width and height numbers in pointings.list to define the
c             FOV footprint. This allows for general polygons, or 40 CCD
c             MegaPrime specific FOV. The routine is still compatible with the
c             old input files.
c
c 2016-05-04: Continue looping on pointings until object is detected,
c             characterized and tracked. Don't stop at first detection.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File generated on 2017-01-12
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      subroutine Detos1 (a, e, inc, node, peri, mt0, jday, hx, color,
     $  gb, ph, period, amp, surnam, seed, flag, ra, dec, d_ra, d_dec,
     $  r, delta, m_int, m_rand, eff, isur, mt, jdayp, ic, surna,
     $  h_rand)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine determines if a given object is seen by the survey
c described in the directory \verb|surnam|.
c An object is described by its ecliptic (J2000) barycentric osculating
c elements given at time \verb|jday|.
c This version uses polygons to describe the footprint of the block on
c the sky.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : January 2006
c Version 2 : May 2013
c Version 3 : March 2016
c Version 4 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file).
c             Continue looping on pointings until object is detected,
c             characterized and tracked. Don't stop at first detection.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     a     : Semi-major axis [AU] (R8)
c     e     : Eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     mt0   : Mean anomaly [rad] (R8)
c     jday  : Time of elements [JD] (R8)
c     hx    : Absolute magnitude of object in 'x' band, what ever this is (R8)
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
c     gb    : opposition surge factor G, Bowell formalism (R8)
c     ph    : phase of lightcurve at epoch jday [rad] (R8)
c     period: period of lightcurve [day] (R8)
c     amp   : amplitude of lightcurve [mag] (R8)
c     surnam: Survey directory name (CH)
c
c OUTPUT
c     seed  : Random number generator seed (I4)
c     flag  : Return flag (I4): 
c                0: not found 
c                1: found, but not tracked
c                2: found and tracked
c                3: characterized, but not tracked
c                4: characterized and tracked
c     ra    : Right ascension at detection [rad] (R8)
c     dec   : Declination at detection [rad] (R8)
c     d_ra  : Right ascension rate [rad/day] (R8)
c     d_dec : Declination rate [rad/day] (R8)
c     r     : Sun-object distance [AU] (R8)
c     delta : Earth-object distance [AU] (R8)
c     m_int : Intrinsic apparent magnitude, in x-band (R8)
c     m_rand: Averaged randomized magnitude, in detection filter (R8)
c     eff   : Actual efficiency of detection (R8)
c     isur  : Identification number of survey the object was in (I4)
c     mt    : Mean anomaly at discovery [rad] (R8)
c     jdayp : Time of discovery [JD] (R8)
c     ic    : Index of color used for survey (I4)
c     surna : Detection survey name (CH10)
c     h_rand: Absolute randomized magnitude, in detection filter (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) a
Cf2py intent(in) e
Cf2py intent(in) inc
Cf2py intent(in) node
Cf2py intent(in) peri
Cf2py intent(in) mt0
Cf2py intent(in) jday
Cf2py intent(in) hx
Cf2py intent(in) color
Cf2py intent(in) gb
Cf2py intent(in) ph
Cf2py intent(in) period
Cf2py intent(in) amp
Cf2py intent(in) surnam
Cf2py intent(in,out) seed
Cf2py intent(out) flag
Cf2py intent(out) ra
Cf2py intent(out) dec
Cf2py intent(out) d_ra
Cf2py intent(out) d_dec
Cf2py intent(out) r
Cf2py intent(out) delta
Cf2py intent(out) m_int
Cf2py intent(out) m_rand
Cf2py intent(out) eff
Cf2py intent(out) isur
Cf2py intent(out) mt
Cf2py intent(out) jdayp
Cf2py intent(out) ic
Cf2py intent(out) surna
Cf2py intent(out) h_rand


      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      integer*4
     $  screen, keybd, verbose, lun_s, lun_h

      parameter
     $  (screen = 6, keybd = 5, verbose = 9,
     $  lun_s = 13, lun_h = 6)

      real*8
     $  a, e, inc, node, peri, mt, hx, h, jday, pos(3), m_int, m_rand,
     $  r, delta, alpha, ra, dec, pos2(3), r2, ra2, dec2, h_rand,
     $  sur_pl(2,n_e_max,n_sur_max),
     $  sur_jday(n_sur_max), sur_ff(n_sur_max),
     $  sur_x(n_sur_max), sur_y(n_sur_max), sur_z(n_sur_max),
     $  sur_r(n_sur_max), sur_x2(n_sur_max), sur_y2(n_sur_max),
     $  sur_z2(n_sur_max), sur_r2(n_sur_max), sur_jday2(n_sur_max),
     $  sur_eff_b(n_bin_max, n_r_max, n_sur_max),
     $  sur_eff_m(n_bin_max, n_r_max, n_sur_max),
     $  sur_mmag(n_sur_max), sur_rt(2, n_r_max, n_sur_max),
     $  sur_rn(n_sur_max), sur_rx(n_sur_max), sur_an(n_sur_max),
     $  sur_aw(n_sur_max), sur_tm(n_sur_max), sur_ts(n_sur_max),
     $  sur_ta(n_sur_max), sur_dm(6,n_sur_max), sur_ph(3,n_sur_max),
     $  sur_ml(n_r_max, n_sur_max), mag_err(6), photf(3), color(*),
     $  maglim, width, height, ra_p, dec_p, jdayp, ff,
     $  obspos(3), ros, obspos2(3), ros2, mag_max, mag_faint, random,
     $  eta, track, jday_o, mt0, jdayp2, tmp, eff_lim,
     $  d_ra, d_dec, r_min, r_max, ang, ang_w,
     $  track_max, track_mag, track_slope, angle, rate, delta2,
     $  mag_peri, dmag, eff, ran3, gb, ph, period, amp,
     $  poly(2,n_e_max), p(2), ra_l, dec_l, d_ra_l, d_dec_l, r_l,
     $  delta_l, m_int_l, m_rand_l, eff_l, mt_l, jdayp_l

      integer*4
     $  sur_code(n_sur_max), sur_eff_n(n_r_max, n_sur_max),
     $  sur_nr(n_sur_max), sur_ne(n_sur_max), sur_f(n_sur_max), i,
     $  filt_i, ic, n_e, point_in_polygon, flag_l,
     $  n_sur, i_sur, ierr, seed, incode, outcod, flag, isur, nph

      character
     $  surnam*(*), sur_eff(n_sur_max)*80, stra*13, stdec*13, surna*10

      logical
     $  debug, newpos, rate_ok, first

      external
     $  eta, ran3

      data
     $  first /.true./,
     $  debug /.false./,
     $  eff_lim /0.4d0/

      save

      flag = 0
      flag_l = 0

      if (first) then
         first = .false.

c Opens and reads in survey definitions
         call GetSurvey (surnam, lun_s,
     $     n_sur, sur_pl, sur_ne, sur_jday, sur_ff,
     $     sur_code, sur_x, sur_y, sur_z, sur_r, sur_jday2,
     $     sur_x2, sur_y2, sur_z2, sur_r2, sur_eff, sur_nr, sur_rt,
     $     sur_eff_n, sur_eff_b, sur_eff_m, sur_mmag, sur_rn, sur_rx,
     $     sur_an, sur_aw, sur_ta, sur_tm, sur_ts, sur_dm, sur_ph,
     $     sur_ml, sur_f, ierr)
         if (ierr .ne. 0) then
            if (ierr .eq. 100) then
               write (screen, *)
     $           'GetSurvey: reached maximum number of pointings, ',
     $           n_sur
            else if (ierr .eq. 10) then
               write (screen, *)
     $           'Unable to open survey file in ',surnam
            else if (ierr .eq. 30) then
               goto 100
            else
               write (screen, *) 'Unknown return code in read_sur.'
            end if
            stop
         end if
 100     continue
c Determine overall faintest 'x' magnitude for all surveys
         mag_faint = 0.d0
         do i_sur = 1, n_sur
c sur_mmag(i_sur) in survey's filter
c mag_max in 'x' filter
            mag_max = sur_mmag(i_sur) - color(sur_f(i_sur))
            if (mag_max .gt. mag_faint) mag_faint = mag_max
            if (debug) then
               write (verbose, *), i_sur, mag_max, mag_faint
            end if
         end do
         if (debug) then
            write (verbose, *) 'Faintest magnitude = ', mag_faint
         end if
c         stop
      end if

c Compute approximate maximum apparent 'x' magnitude
      r_l = a*(1.d0 - e)
      h = hx - amp*0.5d0
c mag_peri in 'x' filter
      call AppMag (r_l, r_l-1.d0, 1.d0, h, gb, alpha, mag_peri, ierr)

      if (mag_peri .le. mag_faint) then
         jday_o = -1.d30

c loop on surveys
         do i_sur = 1, n_sur
            obspos(1) = sur_x(i_sur)
            obspos(2) = sur_y(i_sur)
            obspos(3) = sur_z(i_sur)
            ros = sur_r(i_sur)
            jdayp_l = sur_jday(i_sur)
            obspos2(1) = sur_x2(i_sur)
            obspos2(2) = sur_y2(i_sur)
            obspos2(3) = sur_z2(i_sur)
            ros2 = sur_r2(i_sur)
            jdayp2 = sur_jday2(i_sur)
            ff = sur_ff(i_sur)
            r_min = sur_rn(i_sur)
            r_max = sur_rx(i_sur)
            ang = sur_an(i_sur)
            ang_w = sur_aw(i_sur)
            track_max = sur_ta(i_sur)
            track_mag = sur_tm(i_sur)
            track_slope = sur_ts(i_sur)
            filt_i = sur_f(i_sur)
            do i = 1, 6
               mag_err(i) = sur_dm(i,i_sur)
            end do
            do i = 1, 3
               photf(i) = sur_ph(i,i_sur)
            end do
            n_e = sur_ne(i_sur)
            do i = 1, n_e+1
               poly(1,i) = sur_pl(1,i,i_sur)
               poly(2,i) = sur_pl(2,i,i_sur)
            end do

c Quick and dirty trick to avoid some objects on not too faint surveys:
c drop objects that are fainter, at pericenter, and as seen from the
c Sun, than the faintest magnitude recorded for that survey, in 'x'
c band.
c mag_max in 'x' filter
            mag_max = sur_mmag(i_sur) - color(filt_i)

c Any chance this survey can see the object ?
c mag_peri in 'x' filter
            if (mag_peri .le. mag_max) then

               newpos = .false.
               if (abs(jdayp_l-jday_o) .gt. 0.1d0) then
                  mt_l = mt0
     $              + (twopi/(a**1.5d0*365.25d0))*(jdayp_l-jday)
                  mt_l = mt_l - int(mt_l/twopi)*twopi
                  call pos_cart (a, e, inc, node, peri, mt_l, pos(1),
     $              pos(2), pos(3))
                  jday_o = jdayp_l
                  newpos = .true.
               end if
               if (debug) then
                  write (verbose, *) 'Survey: ', i_sur
                  write (verbose, *) pos(1), pos(2), pos(3), jday,
     $              jday_o
                  write (verbose, *) obspos(1), obspos(2), obspos(3),
     $              jdayp_l
                  write (verbose, *) obspos2(1), obspos2(2), obspos2(3),
     $              jdayp2
               end if
               call DistSunEcl (jdayp_l, pos, r_l)
               call RADECeclXV (pos, obspos, delta_l, ra_l, dec_l)
               p(1) = ra_l
               p(2) = dec_l
c Get mag in actual survey filter.
               h = hx + color(filt_i)
     $           + amp*0.5d0*sin((jdayp_l-jday)/period*twopi+ph)
c mag in survey's filter
               call AppMag (r_l, delta_l, ros, h, gb, alpha, m_int_l,
     $           ierr)
               if (ierr .ne. 0) then
                  write (screen, *) 'AppMag: something''s wrong !'
                  write (screen, *) 'ierr = :', ierr
                  write (screen, *) 'Survey number: ', i_sur
                  stop
               end if

c Format angles for output
               if (debug) then
                  incode = 1
                  outcod = 1
                  call Format (ra_l, incode, outcod, stra, ierr)
                  if (ierr .ne. 0) then
                     write (screen, *) 'Error in formatting output.'
                     write (screen, *) 'ierr = ', ierr
                     stop
                  end if
                  outcod = 0
                  call Format (dec_l, incode, outcod, stdec, ierr)
                  if (ierr .ne. 0) then
                     write (screen, *) 'Error in formatting output.'
                     write (screen, *) 'ierr = ', ierr
                     stop
                  end if
                  write (verbose,
     $              '(3(f8.3, 1x), a13, 1x, a13)')
     $              mt0/drad, peri/drad, node/drad,
     $              stra, stdec
                  write (verbose, *) ra_l/drad, dec_l/drad, m_int_l,
     $              mag_max
               end if

c Still any chance to see it (comparison in survey filter band) ?
c sur_mmag(i_sur) in survey's filter
c mag in survey's filter
               if (m_int_l .le. sur_mmag(i_sur)) then

c Is the object in the FOV ?
c
c Here we use polygons.
                  ierr = point_in_polygon(p, poly, n_e)
                  if (debug) then
                     write (verbose, *) 'Check for FOV.'
                     write (verbose, *) n_e, ierr
                     do i = 1, n_e+1
                        write (verbose, *) poly(1,i)/drad, poly(2,i)
     $                    /drad
                     end do
                  end if
                  if (ierr .gt. 0) then

c Check for chip gaps, ..., the filling factor.
                     random = ran3(seed)
                     if (random .le. ff) then

                        if (debug) then
                           write (verbose, *)
     $                       'In FOV of survey. Check filling factor.'
                           write (verbose, *) random, ff
                        end if
c Well, how is its rate of motion ? Within the rate cut or not ?
                        mt_l = mt0 + (twopi/(a**1.5d0*365.25d0))*(jday_o
     $                    + jdayp2 - jdayp_l - jday)
                        mt_l = mt_l - int(mt_l/twopi)*twopi
                        call pos_cart (a, e, inc, node, peri, mt_l,
     $                    pos2(1), pos2(2), pos2(3))
                        call DistSunEcl (jdayp2, pos2, r2)
                        call RADECeclXV (pos2, obspos2, delta2, ra2,
     $                    dec2)
                        if (debug) then
                           write (verbose, *)
     $                       'Check for second position.'
                           write (verbose, *) mt_l
                           write (verbose, *) pos2(1), pos2(2), pos2(3)
                           write (verbose, *) delta2, ra2, dec2
                        end if
                        d_ra_l = ra_l - ra2
                        if (d_ra_l .gt. Pi) d_ra_l = d_ra_l - TwoPi
                        if (d_ra_l .lt. -Pi) d_ra_l = d_ra_l + TwoPi
                        d_ra_l = d_ra_l/(jdayp2 - jdayp_l)*dcos(dec_l)
                        d_dec_l = (dec2 - dec_l)/(jdayp2 - jdayp_l)
                        rate = dsqrt(d_ra_l**2 + d_dec_l**2)
                        angle = atan2(d_dec_l/rate, d_ra_l/rate)
                        if (angle .lt. -Pi) angle = angle + TwoPi
                        if (angle .gt. Pi) angle = angle - TwoPi
                        rate_ok = (rate .ge. r_min)
     $                    .and. (rate .le. r_max)
                        rate_ok = rate_ok .and.
     $                    (dabs(ang - angle) .le. ang_w)
                        if (debug) then
                           write (verbose, *) 'Check for rate.'
                           write (verbose, *) rate/drad*3600.d0/24.d0,
     $                       r_min/drad*3600.d0/24.d0,
     $                       r_max/drad*3600.d0/24.d0
                           write (verbose, *) angle/drad,
     $                       ang/drad, ang_w/drad
                           write (verbose, *) pos2(1), pos2(2), pos2(3)
                        end if
                        if (rate_ok) then

c Now check for the efficiency
                           eff_l = eta(sur_rt(1,1,i_sur),
     $                       sur_nr(i_sur), sur_eff_n(1,i_sur),
     $                       sur_eff_b(1,1,i_sur), sur_eff_m(1,1,i_sur),
     $                       m_int_l, rate, sur_ml(1,i_sur), maglim)
                           random = ran3(seed)
                           if (debug) then
                              write (verbose, *)
     $                          'Rate OK. Check detection.'
                              write (verbose, *) random, eff_l, maglim
                           end if
                           if (random .le. eff_l) then
c Compute "measured" magnitude with 1 to 3 averaged values
                              random = ran3(seed)
                              call magran (m_int_l, mag_err, seed, tmp,
     $                          dmag)
                              m_rand_l = tmp
                              if (random .gt. photf(1)) then
                                 call magran (m_int_l, mag_err, seed,
     $                             tmp, dmag)
                                 m_rand_l = (m_rand_l + tmp)/2.d0
                              end if
                              if (random .gt. photf(1)+photf(2)) then
                                 call magran (m_int_l, mag_err, seed,
     $                             tmp, dmag)
                                 m_rand_l = (2.d0*m_rand_l + tmp)/3.d0
                              end if
c Determine efficiency of detection for that magnitude
                              eff_l = eta(sur_rt(1,1,i_sur),
     $                          sur_nr(i_sur), sur_eff_n(1,i_sur),
     $                          sur_eff_b(1,1,i_sur),
     $                          sur_eff_m(1,1,i_sur), m_rand_l, rate,
     $                          sur_ml(1,i_sur), maglim)
c Hurray ! We found it.
                              flag_l = 1
                              if (debug) then
                                 write (verbose, *)
     $                             'Hurray ! We found it.'
                              end if

c Determine if tracked
                              random = ran3(seed)
                              track = min(track_max,
     $                          1.d0 + (m_rand_l - track_mag)
     $                          *track_slope)
                              if (debug) then
                                 write (verbose, *)
     $                             'Checking for track. ', random,
     $                             track
                              end if
                              if (random .le. track) then
                                 flag_l = 2
                              end if
c Decide if characterized or not
                              if (maglim .gt. 0.d0) then
                                 if (m_rand_l .le. maglim)
     $                             flag_l = flag_l + 2
                              else
                                 if (eff_l .ge. eff_lim)
     $                             flag_l = flag_l + 2
                              end if
c Record what needs to be recorded.
                              if (flag_l .gt. flag) then
                                 isur = i_sur
                                 surna = sur_eff(i_sur)
     $                             (1:min(len(surna),len(sur_eff(1))))
c Converting intrinsic magnitude to 'x' band, keeping apparent
c magnitude in discovery filter
                                 ic = filt_i
                                 m_int = m_int_l - color(ic)
                                 m_rand = m_rand_l
                                 r = r_l
                                 delta = delta_l
                                 call AbsMag (r, delta, ros, m_rand, gb,
     $                             alpha, h_rand, ierr)
                                 if (ierr .ne. 0) then
                                    write (screen, *)
     $                                'AbsMag: something''s wrong !'
                                    write (screen, *) 'ierr = :', ierr
                                    write (screen, *) 'Survey number: ',
     $                                i_sur
                                    stop
                                 end if
                                 if (debug) then
                                    write (verbose, *)
     $                                'All is good, h_rand.'
                                    write (verbose, *) r, delta
                                    write (verbose, *)
     $                                m_rand, alpha, h_rand
                                 end if
                                 flag = flag_l
                                 ra = ra_l
                                 dec = dec_l
                                 d_ra = d_ra_l
                                 d_dec = d_dec_l
                                 eff = eff_l
                                 mt = mt_l
                                 jdayp = jdayp_l
                              end if
c We got it, and we know if it was tracked and/or characterized.
c Return if tracked and characterized, otherwise keep looping.
                              if (flag .ge. 4) return
                           else
c                              write (6, *) 'Low efficiency: ', a, eff_l,
c     $                          random, m_int_l, rate
                           end if
                        else
c                           write (6, *) 'Rate out of range: ', a, r_min,
c     $                       r_max, rate, ang_w, dabs(ang - angle)
                        end if
                     else
c                        write (6, *) 'Falling in chip gaps: ', a,
c     $                    ff, random
                     end if
                  end if
               else
c                  write (6, *) 'Too faint for this survey: ', a, i_sur,
c     $              sur_mmag(i_sur), m_int_l, hx, filt_i, color(filt_i)
               end if
            else
c               write (6, *) 'Too faint (peri) for this survey: ', a,
c     $           i_sur, mag_max, mag_peri, hx, filt_i, color(filt_i)
            end if

c End loop on surveys
         end do
      else
c         write (6, *) 'Too faint (peri) for all surveys: ', a,
c     $     mag_faint, mag_peri, hx, filt_i, color(filt_i)
      end if

      return

      end

      subroutine AppMag (r, delta, robs, h, g, alpha, mag, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes phase angle and apparent magnitude.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     r     : Sun-object distance (R8)
c     delta : Earth-object distance (R8)
c     robs  : Sun-Earth distance (R8)
c     h     : Absolute magnitude of object (R8)
c     g     : Slope of object (R8)
c
c OUTPUT
c     alpha : Phase angle (R8)
c     mag   : Apparent magnitude (R8)
c     ierr  : Error code
c                0 : nominal run
c               10 : wrong input data
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  ierr

      real*8
     $  r, delta, robs, h, g, alpha, mag, Pi, raddeg, denom, phi1, phi2

      parameter
     $  (Pi = 3.141592653589793d0, raddeg = 180.0d0/Pi)

      ierr = 0
      denom = 2.d0*r*delta
      if (denom .eq. 0.d0) then
         ierr = 10
         return
      end if
      alpha = dacos(dmin1((-robs**2 + delta**2 + r**2)/denom,1.d0))
      phi1 = exp(-3.33d0*(dtan(alpha/2.0d0))**0.63d0)
      phi2 = exp(-1.87d0*(dtan(alpha/2.0d0))**1.22d0)
      mag = 5.d0*dlog10(r*delta) + h
     $  - 2.5d0*dlog10((1.d0 - g)*phi1 + g*phi2)
c      write (6, '(7(f10.4,1x))')
c     $  r, delta, robs, alpha, h, mag,
c     $  2.5d0*dlog10((1.d0 - g)*phi1 + g*phi2)

      return
      end

      subroutine AbsMag (r, delta, robs, mag, g, alpha, h, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes phase angle and absolute magnitude.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2014
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     r     : Sun-object distance (R8)
c     delta : Earth-object distance (R8)
c     robs  : Sun-Earth distance (R8)
c     mag   : Apparent magnitude (R8)
c     g     : Slope of object (R8)
c
c OUTPUT
c     alpha : Phase angle (R8)
c     h     : Absolute magnitude of object (R8)
c     ierr  : Error code
c                0 : nominal run
c               10 : wrong input data
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  ierr

      real*8
     $  r, delta, robs, h, g, alpha, mag, mag0

      h = 0.d0
      call AppMag (r, delta, robs, h, g, alpha, mag0, ierr)
      if (ierr .ne. 0) then
         write (6, *) 'AppMag: something''s wrong !'
         write (6, *) 'ierr = :', ierr
         stop
      end if
      h = mag - mag0

      return
      end

      subroutine Format (angle, incode, outcod, string, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine formats an angle (in rd) into deg, min, sec or hour, min,
c sec. Output is a string.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : September 2003
c Version 2 : March 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     angle : angle to format, radians (r8)
c     incode: 1 input in radian; 0 input in decimal degrees (i)
c     outcod: 1 converts to hours,min,sec; 0 converts to deg.,min,sec (i)
c
c OUTPUT
c     string: Output string (CH)
c     ierr  : Error code
c                0 : nominal run
c               10 : input data code
c               20 : wrong conversion code
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) angle
Cf2py intent(in) incode
Cf2py intent(in) outcod
Cf2py intent(out) string
Cf2py intent(out) ierr

      implicit none

      integer
     $  deg, mn, ierr, incode, outcod, si

      real*8
     $  angle, sec, Pi, raddeg, rm, w

      character
     $  string*13

      parameter
     $  (Pi = 3.141592653589793d0, raddeg = 180.0d0/Pi)

      ierr = 0
      if (incode .eq. 1) then
         w = angle*raddeg
      else if (incode .eq. 0) then
         w = angle
      else
         ierr = 10
         return
      end if
      if (outcod .eq. 1) then
         w = w/15.d0
      else if (outcod .ne. 0) then
         ierr = 20
         return
      end if
      if (w .lt. 0.d0) then
         si = -1
         w = abs(w)
      else
         si = 1
      end if
      deg = int(w)
      rm = (w - deg)*60.d0
      mn = int(rm)
      sec = (rm - mn)*60.d0
      write (string, '(i4.2, 1x, i2.2, 1x, f5.2)') deg, mn, sec
      if (string(9:9) .eq. ' ') string(9:9) = '0'
      if (string(10:10) .eq. ' ') string(10:10) = '0'
      if (si .eq. -1) then
         if (deg .ge. 100) then
            string(1:1) = '-'
         else
            string(2:2) = '-'
         end if
      end if
      deg = si*deg

      return
      end

      subroutine dgauss (i, y)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives a random value with gaussian probability, with 0
c mean and standard deviation 1.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : January 1990
c Version 2 : January 2007
c             Modified to use RAN3 as random number generator rather
c             than PSALUN because it has a much longer periodicity.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     i     : Seed for random number generator (I4)
c
c OUTPUT
c     y     : Random value (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) i
Cf2py intent(out) y

      integer*4
     $  compte, i

      real*8
     $  pi, x1, x2, y, y1, y2, ran3

      external ran3

      data
     1  compte /0/,
     1  pi /3.14159265358979d0/

      save x1, x2, compte, pi

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
      end

      subroutine magran (mag_t, mag_er, seed, mag, magerr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives a randomized value of the magnitude given the
c theoretical magnitude and parameters to compute the uncertainty.
c
c Version 2
c This works for uncertainties given by the measurement on 1 frame only.
c Shouldn't try to combine several frame to estimate the error as this
c mostly account for zeropoint uncertainty and lightcurve.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : January 2006
c Version 2 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     mag_t : Theoretical magnitude of object (R8)
c     mag_er: Magnitude error parameters (6,n*R8)
c     seed  : Seed for random number generator (I4)
c
c OUTPUT
c     mag   : Randomized magnitude (R8)
c     magerr: Magnitude uncertainty (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) mag_t
Cf2py intent(in) mag_er
Cf2py intent(in,out) seed
Cf2py intent(out) mag
Cf2py intent(out) magerr

      implicit none

      real*8
     $  mag_t, mag_er(*), magerr, mag, tmp, mag_th

      integer*4
     $  seed, i

      mag_th = mag_t
c      tmp = log10(mag_er(2)/mag_er(1))/(mag_er(3)-21.d0)
      if (mag_th .le. 21.d0) then
         magerr = mag_er(1)
      else if (mag_th .le. mag_er(3)) then
c         magerr = mag_er(1)*10.d0**(tmp*(mag_th-21.d0))
         magerr = mag_er(1)*10.d0**(mag_er(2)*(mag_th-21.d0))
      else
c         magerr = max(mag_er(1)*10.d0**(tmp*(mag_er(3)-21.d0))
         magerr = max(mag_er(1)*10.d0**(mag_er(2)*(mag_er(3)-21.d0))
     $     - (mag_th - mag_er(3))*mag_er(4), 0.d0)
      end if
      call dgauss(seed, tmp)
      mag = mag_th + magerr*tmp
c      write (19, *) (mag_er(i), i=1,6)
c      write (19, *) mag_th, magerr, tmp, mag
      if (mag_th .gt. mag_er(5)) then
         mag = mag + (mag_th - mag_er(5))*mag_er(6)
      end if
c      write (19, *) mag
c      write (19, '(4(f6.3, 1x), i10)') mag_th, magerr, tmp, mag, seed

      return

      end


      subroutine LatLong (pos, long, lat, r)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine transforms cartesian coordinates into longitude,
c latitute and distance (almost spherical coordinates). If the input
c cordinates are in ICRF, then one obtains the RA and DEC
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : September 2003
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     pos   : Object's cartesian coordinates
c
c OUTPUT
c     long  : Longitude (RA if ICRF)
c     lat   : Latitude (DEC if ICRF)
c     r     : Distance to center
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) pos
Cf2py intent(out) long
Cf2py intent(out) lat
Cf2py intent(out) r

      implicit none

      real*8
     $  pos(3), lat, long, r, Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi)

      r = dsqrt (pos(1)**2 + pos(2)**2 + pos(3)**2)
      long = datan2(pos(2), pos(1))
      if (long .lt. 0.d0) long = long + TwoPi
      lat = asin(pos(3)/r)

      return
      end

      subroutine RADECeclXV (pos, obspos, delta, ra, dec)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the RA and DEC of an object, defined by its
c barycentric ecliptic cartesian coordinates, with respect to an
c observatory, defined by its ICRF cartesian coordinates.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     pos   : Object barycentric ecliptic cartsian coordinates (3*R8)
c     obspos: Observatory ICRF cartsian coordinates (3*R8)
c
c OUTPUT
c     delta : Distance to observatory (R8)
c     ra    : Right Ascension (R8)
c     dec   : Declination (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) pos
Cf2py intent(in) obspos
Cf2py intent(out) delta
Cf2py intent(out) ra
Cf2py intent(out) dec

      implicit none

      real*8
     $  obspos(3), ra, dec, pos(3), opos(3), delta

      integer*4
     $  ierr

c Compute ICRF cartesian coordinates
      call equat_ecl (-1, pos, opos, ierr)
      if (ierr .ne. 0) then
         write (6, *) 'Problem in conversion ecliptic -> equatorial'
      end if

c Compute RA and DEC
      opos(1) = opos(1) - obspos(1)
      opos(2) = opos(2) - obspos(2)
      opos(3) = opos(3) - obspos(3)
      call LatLong (opos, ra, dec, delta)

      return

      end

      real*8 FUNCTION ran3(idum)
Cf2py intent(in,out) idum
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL*8 FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=abs(MSEED-abs(idum))
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
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
      END

      subroutine zero2pi (var)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c This function resets variable 'var' to be between 0 and 2pi
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c B. Gladman  UBC
c Version 1 : January 2007
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT/OUPUT
c     var   : Variable to reset to be between 0 and 2*Pi (R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in,out) var

      implicit none

c Calling arguments
      real*8  var

cSome values better set up as parameters
      real*8
     $  Pi,                     ! Pi
     $  TwoPi                   ! 2*Pi
c Set the values
      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi)

  771 if (var .gt. TwoPi) then
             var = var - TwoPi
             goto 771
      endif
  772 if (var .lt. 0.0d0) then
             var = var + TwoPi
             goto 772
      endif
      return
      end

      subroutine cal2jul (iyyy, mm, dd, jul)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the Julian Day from time given in
c Year, Month, Day (decimal).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2003
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     iyyy  : Year (I4)
c     mm    : Month (I4)
c     dd    : Decimal Day (R8)
c
c OUTPUT
c     mjd   : Julian Day (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in) iyyy
Cf2py intent(in) mm
Cf2py intent(in) dd
Cf2py intent(out) jul

      implicit none

      integer*4 mm, iyyy, id, julday, juld
      real*8 dd, idfrac, jul

      external julday

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
      end

      integer*4 FUNCTION JULDAY(MM,ID,IYYY)

      implicit none

      integer*4 igreg, mm, id, iyyy, jy, jm, ja

      PARAMETER (IGREG=15+31*(10+12*1582))

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
      END

      subroutine ObjAbs (a, e, inc, node, peri, tperi, jday, mag,
     $  code, gb, alpha, h, ra, dec, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine determines the absolute magnitude of an object and its
c sky position (RA, DEC) given its orbital elements (Berstein &
c Kushalani format), measured magnitude and epoch of observation.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : May 2014
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     a     : Semi-major axis [AU] (R8)
c     e     : Eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     tperi : Time of perihelion passage [JD] (R8)
c     jday  : Epoch of observation [JD] (R8)
c     mag   : Apparent magnitude of object (R8)
c     code  : Observatory code (I4)
c              001 : GAIA
c              002 : Geocentric, Mignard's code
c              500 : Geocentric
c     gb    : opposition surge factor, Bowell formalism (R8)
c
c OUTPUT
c     alpha : Phase angle [rad] (R8)
c     h     : Absolute magnitude of object (R8)
c     ra    : Right Ascension (R8)
c     dec   : Declination (R8)
c     ierr  : Error code (I4)
c                0 : nominal run
c               10 : wrong input data
c              100 : date of call earlier than xjdbeg
c              200 : date of call later   than xjdend
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) a
Cf2py intent(in) e
Cf2py intent(in) inc
Cf2py intent(in) node
Cf2py intent(in) peri
Cf2py intent(in) tperi
Cf2py intent(in) jday
Cf2py intent(in) mag
Cf2py intent(in) gb
Cf2py intent(out) alpha
Cf2py intent(out) h
Cf2py intent(out) ra
Cf2py intent(out) dec
Cf2py intent(out) ierr

      implicit none

      real*8
     $  Pi, TwoPi, drad

      integer*4
     $  screen

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0, screen = 6)

      real*8
     $  a, e, inc, node, peri, tperi, mt, h, jday, pos(3),
     $  r, delta, alpha, ra, dec, obpos(3), ros, gb, tmp(3),
     $  mag

      integer*4
     $  ierr, code

      save

      call ObsPos (code, jday, obpos, tmp, ros, ierr)
      if (ierr .ne. 0) then
         write (screen, *)
     $     'Error while computing observatory''s position.'
         write (screen, *) 'ierr = ', ierr
         return
      end if
      mt = (twopi/(a**1.5d0*365.25d0))*(jday-tperi)
      mt = mt - int(mt/twopi)*twopi
      call pos_cart (a, e, inc, node, peri, mt, pos(1),
     $  pos(2), pos(3))
      call DistSunEcl (jday, pos, r)
      call RADECeclXV (pos, obpos, delta, ra, dec)
      call AbsMag (r, delta, ros, mag, gb, alpha, h, ierr)
      if (ierr .ne. 0) then
         write (screen, *) 'AbsMag: something''s wrong !'
         write (screen, *) 'ierr = :', ierr
         return
      end if

      return
      end
      subroutine GetSurvey (survey, lun_s,
     $  n_sur, sur_pl, sur_ne, sur_t, sur_ff, sur_co, sur_x,
     $  sur_y, sur_z, sur_r, sur_t2, sur_x2, sur_y2, sur_z2, sur_r2,
     $  sur_ef, sur_nr, sur_rt, sur_en, sur_eb, sur_em, sur_mm, sur_rn,
     $  sur_rx, sur_an, sur_aw, sur_ta, sur_tm, sur_ts, sur_dm, sur_ph,
     $  sur_ml, sur_f, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine reads in a survey description.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c Version 2 : October 2004
c Version 3 : January 2006
c Version 4 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     survey: Survey directory name (CH)
c     lun_s : Logical unit for file (I4)
c
c OUTPUT
c     n_sur : Number of pointings read (I4)
c     sur_pl: Array of polygons (2,n_e_max,n*R8)
c     sur_ne: Array of nuumber of edges in polygons (n*I4)
c     sur_t : Time of pointing (n*R8)
c     sur_ff: Filling factor (n*R8)
c     sur_co: Observatory code (n*I4)
c     sur_x : ICRF x coordinate of observatory (n*R8)
c     sur_y : ICRF y coordinate of observatory (n*R8)
c     sur_z : ICRF z coordinate of observatory (n*R8)
c     sur_r : Distance from observatory to Sun (n*R8)
c     sur_t2: Time of pointing 2 hours later (n*R8)
c     sur_x2: ICRF x coordinate of observatory 2 hours later (n*R8)
c     sur_y2: ICRF y coordinate of observatory 2 hours later (n*R8)
c     sur_z2: ICRF z coordinate of observatory 2 hours later (n*R8)
c     sur_r2: Distance from observatory to Sun 2 hours later (n*R8)
c     sur_ef: Name of efficiency function file (n*CH80)
c     sur_nr: Number of efficiency functions per pointing (n*I4)
c     sur_rt: Rates limits for efficiency function ([rad/day]) (2,n*R8)
c     sur_en: Number of bins in efficiency function (n_r_max,n*I4)
c     sur_eb: Bin centers for efficiency function (n_bin_max,n_r_max,n*R8)
c     sur_em: Efficiency at bin center (n_bin_max,n_r_max,n*R8)
c     sur_mm: Limiting magnitude for each survey (n*R8)
c     sur_rn: Lower rate cut ([rad/day]) (n*R8)
c     sur_rx: Upper rate cut ([rad/day]) (n*R8)
c     sur_an: Mean direction of motion ([rad]) (n*R8)
c     sur_aw: Half-width of direction cone ([rad]) (n*R8)
c     sur_ta: Maximum tracking fraction (n*R8)
c     sur_tm: Tracking fraction magnitude intercept (n*R8)
c     sur_ts: Tracking fraction magtnidue slope (n*R8)
c     sur_dm: Magnitude error parameters (6,n*R8)
c     sur_ph: Photometric measurments fractions (3,n*R8)
c     sur_ml: Limiting magnitude of survey (n_r_max,n*R8)
c     sur_f : Filter used for this survey (n*I4)
c     ierr  : Error code (I4)
c                0 : nominal run
c              100 : Maximum number of objects reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) survey
Cf2py intent(in) lun_s
Cf2py intent(out) n_sur
Cf2py intent(out) sur_pl
Cf2py intent(out) sur_ne
Cf2py intent(out) sur_t
Cf2py intent(out) sur_ff
Cf2py intent(out) sur_co
Cf2py intent(out) sur_x
Cf2py intent(out) sur_y
Cf2py intent(out) sur_z
Cf2py intent(out) sur_r
Cf2py intent(out) sur_t2
Cf2py intent(out) sur_x2
Cf2py intent(out) sur_y2
Cf2py intent(out) sur_z2
Cf2py intent(out) sur_r2
Cf2py intent(out) sur_ef
Cf2py intent(out) sur_nr
Cf2py intent(out) sur_rt
Cf2py intent(out) sur_en
Cf2py intent(out) sur_eb
Cf2py intent(out) sur_em
Cf2py intent(out) sur_mm
Cf2py intent(out) sur_rn
Cf2py intent(out) sur_rx
Cf2py intent(out) sur_an
Cf2py intent(out) sur_aw
Cf2py intent(out) sur_ta
Cf2py intent(out) sur_tm
Cf2py intent(out) sur_ts
Cf2py intent(out) sur_dm
Cf2py intent(out) sur_ph
Cf2py intent(out) sur_ml
Cf2py intent(out) sur_f
Cf2py intent(out) ierr

      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      real*8
     $  sur_t(n_sur_max), sur_ff(n_sur_max), sur_x(n_sur_max),
     $  sur_y(n_sur_max), sur_z(n_sur_max), sur_r(n_sur_max),
     $  sur_t2(n_sur_max), sur_x2(n_sur_max), sur_y2(n_sur_max),
     $  sur_z2(n_sur_max), sur_r2(n_sur_max),
     $  sur_eb(n_bin_max, n_r_max, n_sur_max),
     $  sur_em(n_bin_max, n_r_max, n_sur_max),
     $  jday_p, ff, obspos(3), ros,
     $  rates(2, n_r_max), sur_rt(2, n_r_max, n_sur_max),
     $  eff_b(n_bin_max, n_r_max), eff_m(n_bin_max, n_r_max),
     $  sur_mm(n_sur_max), sur_rn(n_sur_max), sur_rx(n_sur_max),
     $  sur_an(n_sur_max), sur_ph(3,n_sur_max),
     $  sur_aw(n_sur_max), sur_ta(n_sur_max), sur_tm(n_sur_max),
     $  sur_ts(n_sur_max), sur_dm(6,n_sur_max), 
     $  sur_ml(n_r_max,n_sur_max), eta, obspos2(3),
     $  ros2, maglim(n_r_max),
     $  jday_p2, rate_c(4), d_mag(6), rate, track(3), photf(3), tmp,
     $  poly(2,n_e_max), sur_pl(2,n_e_max,n_sur_max)

      integer*4
     $  sur_co(n_sur_max), sur_en(n_r_max,n_sur_max), sur_nr(n_sur_max),
     $  sur_f(n_sur_max), nr, j,
     $  eff_n(n_r_max), code, i, lun_s, n_sur, ierr, i1, i2, filt_i,
     $  n_e, sur_ne(n_sur_max)

      character
     $  survey*(*), sur_ef(n_sur_max)*(*), eff_name*80

      logical
     $  finished

      external eta

c Open and read in survey definitions
      call read_file_name (survey, i1, i2, finished, len(survey))
      n_sur = 0
 200  continue
         call read_sur (survey(i1:i2), lun_s, poly, n_e,
     $     jday_p, ff, code, obspos, ros, jday_p2,
     $     obspos2, ros2, eff_name, nr, rates, eff_n, eff_b, eff_m,
     $     rate_c, track, d_mag, photf, maglim, filt_i, ierr)

         if (ierr .ne. 0) then
            if (ierr .eq. 10) then
               write (6, *)
     $           'Unable to open ',survey(i1:i2),'/pointings.list'
            else if (ierr .eq. 20) then
               write (6, *)
     $           'Error reading ',survey(i1:i2),'/pointings.list'
               write (6, *) 'Survey number: ', n_sur
               goto 200
            else if (ierr .eq. 30) then
               goto 300
            else
               write (6, *) 'Unknown return code in read_obj.'
            end if
            stop
         end if

         n_sur = n_sur + 1
         sur_t(n_sur) = jday_p
         sur_ff(n_sur) = ff
         sur_co(n_sur) = code
         sur_x(n_sur) = obspos(1)
         sur_y(n_sur) = obspos(2)
         sur_z(n_sur) = obspos(3)
         sur_r(n_sur) = ros
         sur_t2(n_sur) = jday_p2
         sur_x2(n_sur) = obspos2(1)
         sur_y2(n_sur) = obspos2(2)
         sur_z2(n_sur) = obspos2(3)
         sur_r2(n_sur) = ros2
         sur_rn(n_sur) = rate_c(1)
         sur_rx(n_sur) = rate_c(2)
         sur_an(n_sur) = rate_c(3)
         sur_aw(n_sur) = rate_c(4)
         sur_ta(n_sur) = track(1)
         sur_tm(n_sur) = track(2)
         sur_ts(n_sur) = track(3)
         sur_f(n_sur) = filt_i
         do i = 1, 6
            sur_dm(i,n_sur) = d_mag(i)
         end do
         do i = 1, 3
            sur_ph(i,n_sur) = photf(i)
         end do 
         sur_ef(n_sur) = eff_name
         sur_nr(n_sur) = nr
         sur_ne(n_sur) = n_e
         do j = 1, n_e+1
            sur_pl(1,j,n_sur) = poly(1,j)
            sur_pl(2,j,n_sur) = poly(2,j)
         end do
c         write (18, *) 'Survey number: ', n_sur
c         write (18, *) sur_w(n_sur), sur_h(n_sur), sur_ra(n_sur),
c     $     sur_de(n_sur)
c         write (18, *) sur_t(n_sur), sur_ff(n_sur), sur_co(n_sur)
c         write (18, *) sur_x(n_sur), sur_y(n_sur), sur_z(n_sur),
c     $     sur_r(n_sur)
c         write (18, *) sur_ef(n_sur), sur_nr(n_sur)
         sur_mm(n_sur) = 0.d0
         do j = 1, nr
            sur_rt(1,j,n_sur) = rates(1,j)
            sur_rt(2,j,n_sur) = rates(2,j)
            sur_en(j, n_sur) = eff_n(j)
            sur_ml(j, n_sur) = maglim(j)
c            write (18, *) j, sur_rt(1,j,n_sur), sur_rt(2,j,n_sur),
c     $        sur_en(j,n_sur)
            if (eff_n(j) .gt. 0) then
               do i = 1, eff_n(j)
                  sur_eb(i,j,n_sur) = eff_b(i,j)
                  sur_em(i,j,n_sur) = eff_m(i,j)
c                  write (18, *) j,i,sur_eb(i,j,n_sur),sur_em(i,j,n_sur)
               end do
               sur_mm(n_sur) = max(sur_mm(n_sur), eff_b(eff_n(j),j))
c               write (18, *) sur_mm(n_sur)
            else if (eff_n(j) .eq. -1) then
               do i = 1, 3
                  sur_em(i,j,n_sur) = eff_m(i,j)
               end do
            else if ((eff_n(j) .eq. -2) .or. (eff_n(j) .eq. -4)) then
               do i = 1, 4
                  sur_em(i,j,n_sur) = eff_m(i,j)
               end do
            else if (eff_n(j) .eq. -3) then
               do i = 1, 3
                  sur_em(i,j,n_sur) = eff_m(i,j)
               end do
            else
               write (6, *) 'Got efficiency function type ', eff_n(j), j
               write (6, *) 'Should be >0, -1, -2, -3 or -4.'
               stop 'Something is wrong with this. Aborting.'
            end if
            if (eff_n(j) .lt. 0) then
               rate = 0.5d0*(sur_rt(1,j,n_sur) + sur_rt(2,j,n_sur))
c               write (18, *) 'Rate: ', rate
               ff = 40.d0
 250           continue
                  ff = ff - 0.1d0
                  ros = eta(rates, sur_nr(n_sur),
     $              eff_n, eff_b, eff_m, ff, rate, maglim, tmp)
c                  write (18, *) ff, ros
                  if ((ros .eq. 0.d0) .and. (ff .ge. -0.05d0)) goto 250
 260           continue
               sur_mm(n_sur) = max(sur_mm(n_sur), ff+0.1d0)
            end if
c            write (18, *) n_sur, j, sur_mm(n_sur)
         end do
c         write (18, *) n_sur, sur_mm(n_sur)
         goto 200
 300  continue

      ierr = 0
      return
      end

      subroutine read_file_name (base_name, i1, i2, finished, len)

      implicit none

      integer*4
     $  i1, i2, len

      character
     $  base_name*(*)

      logical
     $  finished

      finished = .false.
      i1 = 1
 100  continue
      if ((base_name(i1:i1) .eq. char(0))
     $  .or. (base_name(i1:i1) .eq. char(9))
     $  .or. (base_name(i1:i1) .eq. ' ')) then
         i1 = i1 + 1
         if (i1 .eq. len) then
            finished = .true.
            return
         end if
         goto 100
      end if
 101  continue

      i2 = i1 + 1
 110  continue
      if ((base_name(i2:i2) .ne. char(0))
     $  .and. (base_name(i2:i2) .ne. char(9))
     $  .and. (base_name(i2:i2) .ne. ' ')) then
         if (i2 .eq. len) goto 111
         i2 = i2 + 1
         goto 110
      end if
      i2 = i2 - 1
 111  continue

      return
      end

      subroutine read_sur (dirn, lun_in, poly, n_e, jday,
     $  ff, code, pos, r, jday2, pos2, r2, efnam, nr, rates, eff_n,
     $  eff_b, eff_m, rate_c, track, d_mag, photf, maglim, filt_i, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine opens and reads in the survey description file.
c Angles are returned in radian.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c Version 2 : October 2004
c Version 3 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     dirn  : Name of directory with survey definition (CH)
c     lun_in: File unit (I4)
c
c OUTPUT
c     poly  : Array of points to define a polygon ((2,n)*R8)
c     n_e   : Number of edges in polygon (I4)
c     jday  : Time of pointing (R8)
c     ff    : Filling factor (R8)
c     code  : Observatory code (I4)
c     pos   : Observatory position (3*R8)
c     r     : Observatory distance to Sun (R8)
c     jday2 : Time of pointing 2 hours later (R8)
c     pos2  : Observatory position 2 hours later (3*R8)
c     r2    : Observatory distance to Sun 2 hours later (R8)
c     efnam : Efficiency file name (CH)
c     nr    : Number of efficiency fucntions (I4)
c     rates : Rates limits for efficiency function ([rad/day]) (2,n*R8)
c     eff_n : Number of bins in efficiency (I4)
c     eff_b : efficiency center of bins (n*R8)
c     eff_m : efficiency of corresponding bin (n*R8)
c     rate_c: rate cut parameters ([rad/day] and [rad]) (4*R8)
c     track : tracking fraction parameters (3*R8)
c     d_mag : Magnitude error parameters (6*R8)
c     photf : Photometric measurments fractions (3*R8)
c     maglim: Limiting magnitude of survey (n*R8)
c     filt_i: filter index (I4)
c     ierr  : Error code
c                0 : nominal run
c               10 : unable to open pointing file
c               20 : error reading record
c               30 : end of file reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) dirn
Cf2py intent(in) lun_in
Cf2py intent(out) poly
Cf2py intent(out) n_e
Cf2py intent(out) jday
Cf2py intent(out) ff
Cf2py intent(out) code
Cf2py intent(out) pos
Cf2py intent(out) r
Cf2py intent(out) jday2
Cf2py intent(out) pos2
Cf2py intent(out) r2
Cf2py intent(out) efnam
Cf2py intent(out) nr
Cf2py intent(out) rates
Cf2py intent(out) eff_n
Cf2py intent(out) eff_b
Cf2py intent(out) eff_m
Cf2py intent(out) rate_c
Cf2py intent(out) track
Cf2py intent(out) d_mag
Cf2py intent(out) photf
Cf2py intent(out) maglim
Cf2py intent(out) filt_i
Cf2py intent(out) ierr

      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      integer*4
     $  n_e

      real*8
     $  w, h, ra, dec, ff, pos(3), eff_b(n_bin_max, n_r_max),
     $  eff_m(n_bin_max, n_r_max),
     $  jday, r, vel(3), photf(3), maglim(n_r_max),
     $  pos2(3), r2, jday2, rate_c(4), track(3), d_mag(6),
     $  rates(2,n_r_max), poly(2,n_e_max)

      integer
     $  lun_in, ierr, j, code, eff_n(n_r_max), nw, lw(nw_max), lun_e,
     $  ierr_e, i1, i2, i3, i4, filt_i, nr

      character
     $  line*100, dirn*(*), efnam*(*), word(nw_max)*80, fname*100

      logical
     $  opened, finished

      data opened /.false./

      save opened

      call read_file_name (dirn, i1, i2, finished, len(dirn))
      ierr = 0
      lun_e = lun_in + 1
      if (.not. opened) then
         line(1:i2-i1+1) = dirn(i1:i2)
         line(i2-i1+2:) = '/pointings.list'
         open (unit=lun_in, file=line, status='old', err=1000)
         opened = .true.
      end if
 1500 continue
      do j = 1, len(line)
         line(j:j) = ' '
      end do
      read (lun_in, '(a)', err=2000, end=3000) line
      if (line(1:1) .eq. '#') goto 1500
      call parse (line, nw_max, nw, word, lw)
      if (word(1)(1:4) .eq. 'ears') then
         if (nw .lt. 7) goto 2000
         j = index(word(2), ':')
         if (j .le. 0) then
            read (word(2), *, err=2000) ra
         else
            call hms (word(2), ra)
            ra = ra*15.d0
         end if
         ra = ra*drad
         call hms (word(3), dec)
         dec = dec*drad
         call create_ears(ra, dec, poly, n_e)
         do j = nw, 4, -1
            word(j+1) = word(j)
            lw(j+1) = lw(j)
         end do
         nw = nw + 1
      else if (word(1)(1:4) .eq. 'poly') then
         if (nw .lt. 8) goto 2000
         read (word(2), *, err=2000) n_e
         j = index(word(3), ':')
         if (j .le. 0) then
            read (word(3), *, err=2000) ra
         else
            call hms (word(3), ra)
            ra = ra*15.d0
         end if
         ra = ra*drad
         call hms (word(4), dec)
         dec = dec*drad
         do j = 1, n_e
            read (lun_in, *, err=2000, end=3000) poly(1,j), poly(2,j)
            poly(1,j) = poly(1,j)*drad
            poly(2,j) = poly(2,j)*drad
         end do
         call create_poly(ra, dec, poly, n_e)
c         write (6, *) 'This feature is not implemented yet.'
c         goto 2000
      else
         if (word(1)(1:4) .eq. 'rect') then
            do j = 2, nw
               word(j-1) = word(j)
               lw(j-1) = lw(j)
            end do
            nw = nw - 1
         end if
         if (nw .lt. 8) goto 2000
         read (word(1), *, err=2000) w
         w = w*drad/2.d0
         read (word(2), *, err=2000) h
         h = h*drad/2.d0
         j = index(word(3), ':')
         if (j .le. 0) then
            read (word(3), *, err=2000) ra
         else
            call hms (word(3), ra)
            ra = ra*15.d0
         end if
         ra = ra*drad
         call hms (word(4), dec)
         dec = dec*drad
         call create_rectangle(w, h, ra, dec, poly, n_e)
      end if
      call check_polygon(poly, n_e)
      read (word(5), *, err=2000) jday
      read (word(6), *, err=2000) ff
      read (word(7), *, err=2000) code

c USE OF SLALIB: need to get longitude, latitude and elevation of
c observatory. This is given by the sla_OBS routine. One then needs to
c get the LST (see documentation on EXPLANATION AND EXAMPLES:
c Ephemerides).

      efnam = word(8)
      call read_file_name (efnam, i3, i4, finished, len(efnam))

c Open and read in efficiency function
      fname(1:i2-i1+2) = dirn(i1:i2)//'/'
      fname(i2-i1+3:) = efnam
      call read_eff (fname, lun_e, eff_b, eff_m, eff_n,
     $  rates, nr, rate_c, d_mag, photf, track, maglim, filt_i, ierr_e)

      if (ierr_e .eq. 10) then
         write (6, *) 'Unable to open '//word(8)
         goto 2000
      else if (ierr_e .eq. 0) then
         goto 1610
      else 
         write (6, *) 'Unknown return code in read_sur.'
         stop
      end if
 1610 continue

c Get rid of unused bins at high magnitude for lookup tables
      do i1 = 1, nr
         if (eff_n(i1) .gt. 0) then
            j = eff_n(i1)
 1700       continue
            if (eff_m(j,i1) .le. 0.d0) then
               j = j - 1
               goto 1700
            end if
            eff_n(i1) = amin0(j+1, eff_n(i1))
         end if
      end do

c Computes observatory position at given jday, in ICRF
      call ObsPos (code, jday, pos, vel, r, ierr_e)
      if (ierr_e .ne. 0) then
         write (6, *) 'Error while computing observatory''s position.'
         write (6, *) 'ierr = ', ierr_e
         goto 2000
      end if

c The same, 2 hours later
      jday2 = jday + TwoHours
      call ObsPos (code, jday2, pos2, vel, r2, ierr_e)
      if (ierr_e .ne. 0) then
         write (6, *) 'Error while computing observatory''s position.'
         write (6, *) 'ierr = ', ierr_e
         goto 2000
      end if
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

      subroutine read_eff (filen, lun_in, bin, eff, eff_n,
     $  rates, nrates, rate_c, mag_er, photf, track, maglim, filt_i,
     $  ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine opens and reads in efficiency file.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c Version 2 : April 2013
c             Changed to read new pointings and efficiency file format
c Version 5 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     filen : object element file name
c     lun_in: File unit
c
c OUTPUT
c     bin   : Magnitude bin center (n_bin_max,n_r_max*R8)
c     eff   : Efficiency for that bin (n_bin_max,n_r_max*R8)
c     eff_n : Number of bins in efficiency (n_r_max*I4)
c               -4 : so-called "square" function
c               -3 : piecewise linear function
c               -2 : double hyperbolic tangent
c               -1 : single hyperbolic tangent
c               >0 : number of bins in lookup table
c     rates : Rates limits for efficiency function ([rad/day]) (2,n_r_max*R8)
c     nrates: Number of efficiency fucntions (I4)
c     rate_c: rate cut parameters ([rad/day] and [rad]) (4*R8)
c     mag_er: Magnitude error parameters (6*R8)
c     photf : Photometric measurments fractions (3*R8)
c     track : tracking fraction parameters (3*R8)
c     maglim: Limiting magnitude of survey (n_r_max*R8)
c     filt_i: filter index (I4)
c                1 : g
c                2 : r
c                3 : i
c                4 : z
c                5 : u
c                6 : B
c                7 : V
c                8 : R
c                9 : I
c     ierr  : Error code
c                0 : nominal run
c               10 : unable to open filen
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py inten(in) filen
Cf2py inten(in) lun_in
Cf2py intent(out) bin
Cf2py intent(out) eff
Cf2py intent(out) eff_n
Cf2py intent(out) rates
Cf2py intent(out) nrates
Cf2py intent(out) rate_c
Cf2py intent(out) mag_er
Cf2py intent(out) photf
Cf2py intent(out) track
Cf2py intent(out) maglim
Cf2py intent(out) filt_i
Cf2py intent(out) ierr

      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      real*8
     $  bin(n_bin_max, n_r_max), eff(n_bin_max, n_r_max), rate_c(4),
     $  mag_er(6), track(3), rates(2, n_r_max), photf(3),
     $  maglim(n_r_max)

      integer
     $  lun_in, ierr, eff_n(n_r_max), eq_ind, nw, lw(nw_max), i, nrates,
     $  j, filt_i

      character
     $  line*100, filen*(*), word(nw_max)*80

      logical
     $  rcut, tr, fi, mag, in_rates, in_func, rate(0:n_r_max), ph

      rcut = .false.
      tr = .false.
      fi = .false.
      mag = .false.
      in_rates = .false.
      in_func = .false.
      rate(0) = .true.
      do i = 1, 10
         rate(i) = .false.
      end do
      ph = .false.

      ierr = 0
      open (unit=lun_in, file=filen, status='old', err=1000)
      nrates = 0

 1500 continue
      read (lun_in, '(a)', err=1500, end=3000) line
      if (line(1:1) .eq. '#') goto 1500
      eq_ind = index(line, '=')
      if (eq_ind .le. 0) goto 1500
      call parse (line(1:eq_ind-1), nw_max, nw, word, lw)
      if (nw .ne. 1) goto 1500
      if (word(1)(1:lw(1)) .eq. 'rate_cut') then
         read (line(eq_ind+1:), *, err=1500) (rate_c(i),i=1,4)
c Change rates to rad/day
         rate_c(1) = rate_c(1)*24.d0/3600.d0*drad
         rate_c(2) = rate_c(2)*24.d0/3600.d0*drad
c Change angles to radian
         rate_c(3) = rate_c(3)*drad
         rate_c(4) = rate_c(4)*drad
         rcut = .true.
         in_rates = .false.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'mag_error') then
         read (line(eq_ind+1:), *, err=1500, end=1500) (mag_er(i),i=1,6)
         mag_er(2) = log10(mag_er(2)/mag_er(1))/(mag_er(3)-21.d0)
         mag = .true.
         in_rates = .false.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'phot_frac') then
         read (line(eq_ind+1:), *, err=1500, end=1500) (photf(i),i=1,3)
         ph = .true.
         in_rates = .false.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'track_frac') then
         read (line(eq_ind+1:), *, err=1500) (track(i),i=1,3)
         tr = .true.
         in_rates = .false.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'filter') then
         call parse (line(eq_ind+1:), nw_max-1, nw, word(2), lw(2))
         if (word(2)(1:1) .eq. 'g') then
            filt_i = 1
         else if (word(2)(1:1) .eq. 'r') then
            filt_i = 2
         else if (word(2)(1:1) .eq. 'i') then
            filt_i = 3
         else if (word(2)(1:1) .eq. 'z') then
            filt_i = 4
         else if (word(2)(1:1) .eq. 'u') then
            filt_i = 5
         else if (word(2)(1:1) .eq. 'B') then
            filt_i = 6
         else if (word(2)(1:1) .eq. 'V') then
            filt_i = 7
         else if (word(2)(1:1) .eq. 'R') then
            filt_i = 8
         else if (word(2)(1:1) .eq. 'I') then
            filt_i = 9
         else
            goto 1500
         end if
         fi = .true.
         in_rates = .false.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'rates') then
         j = nrates
         if (rate(nrates)) then
            j = j + 1
         end if
         read (line(eq_ind+1:), *, err=1500) (rates(i, j),i=1,2)
c Change rates to rad/day
         rates(1,j) = rates(1,j)*24.d0/3600.d0*drad
         rates(2,j) = rates(2,j)*24.d0/3600.d0*drad
         nrates = j
         maglim(nrates) = -1.d0
         in_rates = .true.
         in_func = .false.
      else if (word(1)(1:lw(1)) .eq. 'function') then
         if (in_rates) then
            call parse (line(eq_ind+1:), nw_max-1, nw, word(2), lw(2))
            if (word(2) .eq. 'single') eff_n(nrates) = -1
            if (word(2) .eq. 'double') eff_n(nrates) = -2
            if (word(2) .eq. 'linear') eff_n(nrates) = -3
            if (word(2) .eq. 'square') eff_n(nrates) = -4
            if (word(2) .eq. 'lookup') eff_n(nrates) = 0
            in_func = .true.
         end if
      else if (word(1)(1:lw(1)) .eq. 'linear_param') then
         if (in_func .and. (eff_n(nrates) .eq. -3)) then
            read (line(eq_ind+1:), *, err=1500, end=1500)
     $        (eff(i,nrates),i=1,3)
            rate(nrates) = .true.
c            in_rates = .false.
            in_func = .false.
         end if
      else if (word(1)(1:lw(1)) .eq. 'single_param') then
         if (in_func .and. (eff_n(nrates) .eq. -1)) then
            read (line(eq_ind+1:), *, err=1500, end=1500)
     $        (eff(i,nrates),i=1,3)
            rate(nrates) = .true.
c            in_rates = .false.
            in_func = .false.
         end if
      else if (word(1)(1:lw(1)) .eq. 'double_param') then
         if (in_func .and. (eff_n(nrates) .eq. -2)) then
            read (line(eq_ind+1:), *, err=1500, end=1500)
     $        (eff(i,nrates),i=1,4)
            rate(nrates) = .true.
c            in_rates = .false.
            in_func = .false.
         end if
      else if (word(1)(1:lw(1)) .eq. 'square_param') then
         if (in_func .and. (eff_n(nrates) .eq. -4)) then
            read (line(eq_ind+1:), *, err=1500, end=1500)
     $        (eff(i,nrates),i=1,4)
            rate(nrates) = .true.
c            in_rates = .false.
            in_func = .false.
         end if
      else if (word(1)(1:lw(1)) .eq. 'lookup_param') then
         if (in_func .and. (eff_n(nrates) .ge. 0)) then
            i = eff_n(nrates) + 1
            read (line(eq_ind+1:), *, err=1500, end=1500)
     $        bin(i, nrates), eff(i, nrates)
            eff_n(nrates) = i
            rate(nrates) = .true.
         end if
      else if (word(1)(1:lw(1)) .eq. 'mag_lim') then
         if (in_rates) then
            read (line(eq_ind+1:), *, err=1500, end=1500) maglim(nrates)
         end if
      else
         write (6, *) 'WARNING: unknown key '//word(1)(1:lw(1))
         in_rates = .false.
         in_func = .false.
      end if
      goto 1500

 1000 continue
      ierr = 10
      return

 3000 continue
      close (lun_in)
      if (.not. rcut) then
         write (6, *) 'Survey file: ', filen
         write (6, *) 'ERROR: rate cut parameters not defined.'
         stop
      end if
      if (.not. tr) then
         write (6, *) 'Survey file: ', filen
         write (6, *) 'ERROR: tracking parameters not defined.'
         stop
      end if
      if (.not. fi) then
         write (6, *) 'Survey file: ', filen
         write (6, *) 'ERROR: filter not defined.'
         stop
      end if
 3010 continue
      if ((nrates .gt. 0) .and. .not. rate(nrates)) then
         nrates = nrates - 1
         goto 3010
      end if
      if (nrates .le. 0) then
         write (6, *) 'Survey file: ', filen
         write (6, *) 'ERROR: no efficiency function defined.'
         stop
      end if
      if (.not. mag) then
         mag_er(1) = 0.026d0
         mag_er(2) = 0.33d0
         mag_er(3) = 24.45d0
         mag_er(4) = 0.7d0
         mag_er(5) = 23.7d0
         mag_er(6) = -0.3d0
         write (6, *) 'Survey file: ', filen
         write (6, *)
     $     'WARNING: magnitude error parameters not defined, using '
     $     //'default values:', mag_er(1), ', ', mag_er(2), ', ',
     $     mag_er(3), ', ', mag_er(4), ', ',  mag_er(5), ', ', mag_er(6)
         mag_er(2) = log10(mag_er(2)/mag_er(1))/(mag_er(3)-21.d0)
      end if
      if (.not. ph) then
         photf(1) = 1.d0
         photf(2) = 0.d0
         photf(3) = 0.d0
         write (6, *) 'Survey file: ', filen
         write (6, *)
     $     'WARNING: photometric measurements fractions not defined, '
     $     //'using default values:', photf(1), ', ', photf(2),
     $     ', ', photf(3)
      end if
      do j = 1, nrates
         if (rate(j) .and. (maglim(j) .le. 0.d0)) then
            write (6, *) 'Survey file: ', filen, '; rate range: ',
     $        rates(1,j), ' - ', rates(2,j)
            write (6, *)
     $        'WARNING: limiting magnitude not defined, '
     $        //'using efficiency instead.'
         end if
      end do
      ierr = 0
      return

      end

      subroutine create_ears(ra, dec, poly, n_e)

      implicit none

      integer*4
     $  n_e

      real*8
     $  ra, dec, dra, h, w, poly(2,*), e_w, e_h

c Below values are from assuming a size of 1°x1° for the "central"
c square and add 1/2° height and 1/9° width ears. However, looking at
c the resulting polygons, they seem to be somewhat too large
c      data
c     $  h /0.008726646259971648d0/,
c     $  w /0.008726646259971648d0/,
c     $  e_h /0.004363323129985824d0/,
c     $  e_w /0.0019392547244381439d0/
c
c Now, maybe we can do better. Let's look at image 1805373p.fits.
c Pointing is (13:32:29.29; -9:31:03.9) or (203.122042; -9.517750)
c Now, the corners of the "central" square are:
c (13:30:29.97; -10:00:48.4) or (202.624875; -10.013444)
c (13:34:27.66; -10:00:36.0) or (203.615250; -10.010000)
c (13:34:26.45; -09:00:58.4) or (203.610208; -09.016222)
c (13:30:29.44; -09:01:11.4) or (202.622667; -09.019833)
c This is a "square" of 0.9753° x 0.9937°
c
c Now, add the ears:
c (13:34:53.45; -09:16:26.2) or (203.722708; -09.273944)
c (13:34:27.30; -09:16:26.1) or (203.613750; -09.273917)
c (13:34:27.89; -09:45:09.2) or (203.616208; -09.752556)
c (13:34:54.05; -09:45:04.4) or (203.725208; -09.751222)
c This is a "rectangle" of 0.1075° x 0.4780°
c
c (13:30:29.09; -09:16:38.9) or (202.621208; -09.277472)
c (13:30:02.88; -09:16:42.3) or (202.512000; -09.278417)
c (13:30:03.18; -09:45:19.7) or (202.513250; -09.755472)
c (13:30:29.33; -09:45:21.5) or (202.622208; -09.755972)
c This is a "rectangle" of 0.1076° x 0.4778°
c
c So let's make the area 0.9753° x 0.9937° + 2 x 0.010755° x 0.4779° =
c 0.979435239 sq.deg. The CCDs are 2048x4612 with pixels
c 0.18689x0.18689" resulting in an area of 0.0254558 sq.deg. for each
c CCD, or 1.0182 sq.deg. Oops, there is a big problem. I guess it's the
c size of the pixels that's to big. Applying this size to the cenral
c square yield something reasonable, with pixels covering less than the
c whole size of the array. But the ears are in trouble. The pixel heigh
c seems to be too large. Actually, the main problem is that the WCS is
c not good enough, and the central horizontal gap has vanished, and
c some pixels are enven overlapping. I need to get an image with
c Stephen's header.
c Officially, the header says the size of the pixel is 0.185"x0.185",
c somewhat smaller than what JJ says, but I'll stick to what JJ says.
C Actually, I cannot use Stephen's header with DS9 as the latter
C doesn't know how to use the PVs, and is rather using the CDs.
c
c The horizontal size of the central square shuold be greater than
c 9x2112 = 19008 pixels or 3552.4". This implies a gap of ~32 pixels
c between the chips, in addition to the oversans, 32 pixels on each
c side. I'll assume the small horizontal gap is similar in the center
c of the frame. Then for the same 1° full size, the wisth of the large
c gaps is 327 and 328 pixels.
c
c From all this, I assume a half width of 0.5°, a half height of 0.5°,
c the width of the ears (32+2112)*0.18689" = 1/9°, and half height of
c (2*4644+32)*0.18689" = 0.483837°
c      data
c     $  h /0.008726646259971648d0/,
c     $  w /0.008726646259971648d0/,
c     $  e_h /0.004222278224995352d0/,
c     $  e_w /0.0019392547244381439d0/
c
c The above gives a lot of overlap. Using the PVs and CDs, I determined
c the exact footprint of a series of MegaPrime40 frames (in ~/Research
c /OSSOS/src/MegaPrime40.FOV) and averaged them. In addition to the
c positions below, I also determined the pixel size: 0.186" x 0.1849"
      data
     $  h /0.008678d0/,
     $  w /0.008545d0/,
     $  e_h /0.004169d0/,
     $  e_w /0.001913d0/

      n_e = 12
      poly(2,1) = dec - h
      poly(2,12) = poly(2,1)
      poly(2,13) = poly(2,1)
      poly(2,6) = dec + h
      poly(2,7) = poly(2,6)
      poly(2,2) = dec - e_h
      poly(2,3) = poly(2,2)
      poly(2,10) = poly(2,2)
      poly(2,11) = poly(2,2)
      poly(2,4) = dec + e_h
      poly(2,5) = poly(2,4)
      poly(2,8) = poly(2,4)
      poly(2,9) = poly(2,4)

      dra = w/dcos(poly(2,1))
      poly(1,1) = ra - dra
      poly(1,13) = poly(1,1)
      poly(1,12) = ra + dra
      dra = w/dcos(poly(2,2))
      poly(1,2) = ra - dra
      poly(1,11) = ra + dra
      dra = (w + e_w)/dcos(poly(2,2))
      poly(1,3) = ra - dra
      poly(1,10) = ra + dra
      dra = (w + e_w)/dcos(poly(2,4))
      poly(1,4) = ra - dra
      poly(1,9) = ra + dra
      dra = w/dcos(poly(2,4))
      poly(1,5) = ra - dra
      poly(1,8) = ra + dra
      dra = w/dcos(poly(2,6))
      poly(1,6) = ra - dra
      poly(1,7) = ra + dra
      return
      end

      subroutine create_rectangle(w, h, ra, dec, poly, n_e)

      implicit none

      integer*4
     $  n_e

      real*8
     $  w, h, ra, dec, dra, poly(2,*)

      n_e = 4
      poly(2,1) = dec - h
      poly(2,4) = poly(2,1)
      poly(2,2) = dec + h
      poly(2,3) = poly(2,2)
      poly(2,5) = poly(2,1)

      dra = w/dcos(poly(2,1))
      poly(1,1) = ra - dra
      poly(1,4) = ra + dra
      dra = w/dcos(poly(2,2))
      poly(1,2) = ra - dra
      poly(1,3) = ra + dra
      poly(1,5) = poly(1,1)
      return
      end

      subroutine create_poly(ra, dec, poly, n_e)

      implicit none

      integer*4
     $  n_e, j

      real*8
     $  ra, dec, dra, poly(2,*)

      do j = 1, n_e
         poly(2,j) = dec + poly(2,j)
         poly(1,j) = ra + poly(1,j)/dcos(poly(2,j))
      end do
      poly(1,n_e+1) = poly(1,1)
      poly(2,n_e+1) = poly(2,1)

      return
      end

c \subroutine{parse}

c Parses a line returns a list of words by getting rid of space characters

      subroutine parse (command, nwmax, nw, word, lw)

c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|command| = command line to parse \\
c \verb|lw()| = word lengthes \\
c \verb|nw| = number of words in command line \\
c \verb|nwmax| = maximum allowed number of words \\
c \verb|word()| = words in command line
c \end{verse}

c \subsubsection{Declarations}

      integer*4
     $     lw(1), nw, nwmax

      character
     $     command*(*), word(1)*(*)

c \subsection{Variables}
c \subsubsection{Internal variables}
c \begin{verse}
c \verb|k| = dummy index \\
c \verb|lc| = length of command line \\
c \end{verse}

c \subsubsection{Intrinsic Fortran functions used}
c \begin{verse}
c \verb|index|
c \end{verse}

c \subsubsection{Declarations}

      integer*4
     $     k, lc, lw0

c \subsection{Parsing}

      do nw = 1, nwmax
         lw(nw) = 0
      end do
      lc = len(command)
 1000 continue
      if ((command(lc:lc) .eq. char(0))
     $  .or. (command(lc:lc) .eq. ' ')) then
         lc = lc - 1
         if (lc .eq. 0) goto 1001
         goto 1000
      end if
 1001 continue
      nw = 0
      do k = 1, nwmax
          word(k) = ' '
      end do

 1100 continue
      if (lc .gt. 0) then

c Gets rid of leading space characters
         if (nw .ge. nwmax) then
            write (6, *) command
            write (6, *) 'parse: too many words in command line.'
            stop
         end if
 1050    continue
         if (command(1:1) .eq. ' ') then
            command = command (2:lc)
            lc = lc - 1
            goto 1050
         end if

c Finds a word

         nw = nw + 1
         lw0 = index(command, ' ') - 1
         if (lw0 .le. 0) then
            lw(nw) = lc
            word(nw) = command(1:lc)
            lc = -1
         else
            word (nw) = command (1:lw0)
            lw(nw) = lw0
            command = command (lw0+2:lc)
            lc = lc - lw0 - 1
         end if
         goto 1100
      end if

      return
      end

      subroutine  hms(str,val)
c
c...Crack String And Create Value
c
      IMPLICIT
     *            NONE
      CHARACTER*(*)
     *            str
      DOUBLE PRECISION
     *            val, piece(3), dp, sgn, z
      INTEGER
     *            nstr, i, j, dpfind
      CHARACTER*1
     *            c
c
c...Initialization
c
  100 val = 0.0D00
      DO i=1,3
        piece(i) = 0.0D00
      ENDDO
      j = 1
      dpfind = 0
      sgn = 1.0D00
      nstr = LEN(str)
      IF (nstr.le.0) RETURN
c
c...Loop Over The String
c
      DO i=1,nstr
        c = str(i:i)
c
c...Parse
c
        IF ((c.eq.'-').or.(c.eq.'e').or.(c.eq.'E')
     *  .or.(c.eq.'s').or.(c.eq.'S')) THEN
          sgn = -1.0D00
        ELSEIF ((c.eq.'+').or.(c.eq.'w').or.(c.eq.'W')
     *      .or.(c.eq.'n').or.(c.eq.'N')) THEN
          sgn = 1.0D00
        ELSEIF ((c.eq.':').or.(c.eq.',').or.(c.eq.' ')) THEN
          j = j+1
          dpfind = 0
          IF (j.gt.3) GO TO 110
        ELSEIF (c.eq.'.') THEN
          dpfind = 1
          dp = 1.0D00
        ELSEIF ((c.ge.'0').and.(c.le.'9')) THEN
          z = ICHAR(c)-ICHAR('0')
          IF (dpfind.eq.0) THEN
            piece(j) = 10.0D00*piece(j) + z
          ELSE
            dp = 0.1D00*dp
            piece(j) = piece(j) + dp*z
          ENDIF
        ENDIF
      ENDDO
c
c...Return
c
  110 val = piece(1) + piece(2)/60.0D00 + piece(3)/3600.0D00
      val = val*sgn
      RETURN
      END

      real*8 function eta_raw (rates, nr, eff_n, eff_b,
     $  eff_m, mdum, rdum, ml, maglim)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the efficiency at a given magnitude.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : April 2004
c Version 2 : July 2004
c             Added 3 different types of efficiency functions, namely a
c             single and a double hyperbolic tangent, and a piecewise
c             linear function.
c Version 3 : June 2013
c             Added efficiency dependance on rate
c Version 4 : April 2014
c             Added limiting magnitude determination function of rate
c Version 5 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     rates : Rates limits for efficiency function (2,n*R8)
c     nr    : Number of efficiency fucntions (I4)
c     eff_n : Number of efficiency bins (n*I4)
c     eff_b : Magnitude bin center (nb_max,n*R8)
c     eff_m : Efficiency for that magnitude (nb_max,n*R8)
c     mdum  : magnitude (R8)
c     rdum  : rate of motion (R8)
c     ml    : Limiting magnitudes for rate ranges (n*R8)
c
c OUTPUT
c     maglim: Limiting magnitude at given rate (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) rates
Cf2py intent(in) nr
Cf2py intent(in) eff_n
Cf2py intent(in) eff_b
Cf2py intent(in) eff_m
Cf2py intent(in) mdum
Cf2py intent(in) rdum
Cf2py intent(in) ml
Cf2py intent(out) maglim

      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      integer*4
     $  ilo, ihi, i, eff_n(*), nr, ir

      real*8
     $  m, x, tanh, eff_b(n_bin_max,n_r_max), eff_m(n_bin_max,n_r_max),
     $  mdum, rdum, rates(2,n_r_max), r, ml(n_r_max), maglim

      tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))

c Retrieve magnitude and rate of motion
      m = dmax1(0.d0, mdum)
      r = rdum

c Determine which efficiency function to use
      ir = 0
 500  continue
         ir = ir + 1
         if ((r-rates(1,ir))*(r-rates(2,ir)) .le. 0.d0) goto 600
         if (ir .lt. nr) goto 500
 510  continue
      eta_raw = 0.d0
      maglim = 0.d0
      return
 600  continue
c Report limiting magnitude for this rate range
      maglim = ml(ir)

c This is the direct piecewise linear function from a lokup table
      if (eff_n(ir) .gt. 0) then

c If off bins, then flat continuation.
         if (m .lt. eff_b(1,ir)) then
            eta_raw = eff_m(1,ir)
            return
         else if (m .gt. eff_b(eff_n(ir),ir)) then
            eta_raw = 0.d0
            return
         else

c Linear interpolation of table.
            ilo = 1
            ihi = eff_n(ir)
 1000       continue
            if (ihi - ilo .gt. 1) then
               i = (ihi + ilo)/2
               if (eff_b(i,ir) .lt. m) then
                  ilo = i
               else if (eff_b(i,ir) .gt. m) then
                  ihi = i
               else
                  eta_raw = eff_m(i,ir)
                  return
               end if
               goto 1000
            end if
            eta_raw = eff_m(ilo,ir) + (eff_m(ihi,ir) - eff_m(ilo,ir))*
     $        (m - eff_b(ilo,ir))/(eff_b(ihi,ir) - eff_b(ilo,ir))
         end if

c This is a single hyperbolic tangent function.
c \begin{equation}
c (A/2) * (1. - tanh((R-R_c)/d))
c \end{equation}
      else if (eff_n(ir) .eq. -1) then
         eta_raw = eff_m(1,ir)/2.d0 *
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(3,ir)))

c This is a double hyperbolic tangent function.
c \begin{equation}
c (A/4) * (1. - tanh((R-R_c)/d1)) * (1. - tanh((R-R_c)/d2))
c \end{equation}
      else if (eff_n(ir) .eq. -2) then
         eta_raw = eff_m(1,ir)/4.d0 *
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(3,ir))) * 
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(4,ir)))

c This is a piecewize linear function.
c \begin{eqnarray}
c A & {\rm if} & m < R_1 \\
c \frac{(m - R_2) A}{R_1 - R_2} & {\rm if} & R_1 \le m < R_2 \\
c 0 & {\rm if} & m \ge R_2
c \end{eqnarray}
      else if (eff_n(ir) .eq. -3) then
         if (m .lt. eff_m(2,ir)) then
            eta_raw = eff_m(1,ir)
         else if (m .lt. eff_m(3,ir)) then
            eta_raw = (m - eff_m(3,ir))*eff_m(1,ir)/
     $        (eff_m(2,ir) - eff_m(3,ir))
         else
            eta_raw = 0.d0
         end if

c This is the SKADS defined function
c \begin{equation}
c (A - c * (R-21.0)^2) / (1. + exph((R-R_c)/d))
c \end{equation}
      else if (eff_n(ir) .eq. -4) then
         if (m .lt. 21.d0) then
            eta_raw = eff_m(1,ir)
         else
            eta_raw = (eff_m(1,ir) -
     $        eff_m(2,ir) * (m - 21.d0)**2) /
     $        (1.d0 + exp((m - eff_m(3,ir))/eff_m(4,ir)))
         end if

c Unsupported efficiency function type.
      else
         write (6, *) 'Got efficiency function type ', eff_n(ir)
         write (6, *) 'Should be >0, -1, -2 or -3.'
         stop 'Something is wrong with this. Aborting.'
      end if

      return
      end

      real*8 function eta (rates, nr, eff_n, eff_b,
     $  eff_m, mdum, rdum, ml, maglim)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the efficiency at a given magnitude, keeping
c only the part that can be trusted for actual detectability of the
c theoretical magnitude ($\eta > 0.4$).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c Version 3 : June 2013
c             Added efficiency dependance on rate
c Version 4 : April 2014
c             Added limiting magnitude determination function of rate
c Version 5 : May 2016
c             Changed API to remove size of arrays, added parameter
c             statement to define array sizes (in include file)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     rates : Rates limits for efficiency function (2,n*R8)
c     nr    : Number of efficiency fucntions (I4)
c     eff_n : Number of efficiency bins (n*I4)
c     eff_b : Magnitude bin center (nb_max,n*R8)
c     eff_m : Efficiency for that magnitude (nb_max,n*R8)
c     mdum  : magnitude (R8)
c     rdum  : rate of motion (R8)
c     ml    : Limiting magnitudes for rate ranges (n*R8)
c
c OUTPUT
c     maglim: Limiting magnitude at given rate (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) rates
Cf2py intent(in) nr
Cf2py intent(in) eff_n
Cf2py intent(in) eff_b
Cf2py intent(in) eff_m
Cf2py intent(in) mdum
Cf2py intent(in) rdum
Cf2py intent(in) ml
Cf2py intent(out) maglim

      implicit none


      integer*4
     $  n_sur_max, n_bin_max, n_r_max, n_e_max, nw_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  nw_max = 10)

      real*8
     $  Pi, drad, TwoHours, eps, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0,
     $  TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14)

      integer*4
     $  eff_n(n_r_max), nr

      real*8
     $  eta_raw, eff_b(n_bin_max,n_r_max), eff_m(n_bin_max,n_r_max),
     $  mdum, rdum, lim, rates(2,n_r_max), ml(n_r_max), maglim

      external
     $  eta_raw

      data
     $  lim /0.01d0/

      save lim

      eta = eta_raw(rates, nr, eff_n, eff_b, eff_m, mdum, rdum,
     $  ml, maglim)
      if (eta .lt. lim) eta = 0.d0

      return
      end
      subroutine BaryXV (jday, pos, vel, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives the position and velocity in ecliptic heliocentric
c reference frame of the Solar System barycenter at a given time. Uses
c PlanetXV.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     jday  : Time of elements (Julian day) (R8)
c
c OUTPUT
c     pos   : Position vector (3*R8)
c     vel   : Velocity vector (3*R8)
c     ierr  : Error code (I4)
c                0 : nominal run
c               20 : jday out of range
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) jday
Cf2py intent(out) pos
Cf2py intent(out) vel
Cf2py intent(out) ierr

      implicit none

      real*8
     $  jday, pos(3), vel(3), Pi, TwoPi, jday_min,
     $  jday_max, drad, mu

      integer*4
     $  ind, ierr, n_planets, istat

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  jday_min = 2415020.0, jday_max = 2488070.0, n_planets = 8,
     $  mu = TwoPi**2)

      real*8
     $  pos_p(3), vel_p(3), masses(n_planets), msys, mp

      data masses /6023600.0d0,    ! Mercury
     $              408523.7d0,    ! Venus
     $              328900.56d0,   ! Earth + Moon
     $             3098708.0d0,    ! Mars
     $                1047.349d0,  ! Jupiter
     $                3497.898d0,  ! Saturn
     $               22902.98d0,   ! Uranus
     $               19412.24d0/   ! Neptune

      ierr = 0

c Jday out of range.
      if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
         ierr = 20
         return
      end if

c Ok, do the math.
      pos(1) = 0.d0
      pos(2) = 0.d0
      pos(3) = 0.d0
      vel(1) = 0.d0
      vel(2) = 0.d0
      vel(3) = 0.d0
      msys = mu
      do ind = 1, n_planets
         call PlanetXV (ind, jday, pos_p, vel_p, istat)
         if (istat .ne. 0) then
            ierr = istat
            return
         end if
         mp = mu/masses(ind)
         msys = msys + mp
         pos(1) = pos(1) + mp*pos_p(1)
         pos(2) = pos(2) + mp*pos_p(2)
         pos(3) = pos(3) + mp*pos_p(3)
         vel(1) = vel(1) + mp*vel_p(1)
         vel(2) = vel(2) + mp*vel_p(2)
         vel(3) = vel(3) + mp*vel_p(3)
      end do
      pos(1) = pos(1)/msys
      pos(2) = pos(2)/msys
      pos(3) = pos(3)/msys
      vel(1) = vel(1)/msys
      vel(2) = vel(2)/msys
      vel(3) = vel(3)/msys

      return

      end

      subroutine DistSunEcl (jday, pos, r)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes distance from object defined by barycentric
C ecliptic coordinates to Sun at given jday.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     jday  : Time (R8)
c     pos   : Object barycentric ecliptic cartsian coordinates (3*R8)
c
c OUTPUT
c     r     : Distance from object to Sun (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) jday
Cf2py intent(in) pos
Cf2py intent(out) r

      implicit none

      real*8
     $  jday, pos(3), pos_b(3), vel_b(3), r, jday_b

      integer*4
     $  ierr

      logical
     $  done

      data
     $  done /.false./

      save done, pos_b, vel_b, jday_b

      if ((.not. done) .or. (jday .ne. jday_b)) then
         call BaryXV (jday, pos_b, vel_b, ierr)
         if (ierr .ne. 0) then
            write (6, *) 'Problem while getting barycenter position'
         end if
         done = .true.
         jday_b = jday
      end if

      r = dsqrt((pos(1) + pos_b(1))**2 + (pos(2) + pos_b(2))**2
     $  + (pos(3) + pos_b(3))**2)
      return

      end

      subroutine ObsPos (code, t, pos, vel, r, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the cartesian coordinates of the observatory.
c Reference frame : ICRF
c Units : AU
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : September 2003
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     code  : Observatory code (I4)
c              001 : GAIA
c              002 : Geocentric, Mignard's code
c              500 : Geocentric
c     t     : Time of observation (Julian day, not MJD) (R8)
c
c OUTPUT
c     pos   : Cartesian coordinates of observatory (AU) *(R8)
c     vel   : Cartesian velocity of observatory (AU) *(R8)
c     r     : Distance from observatory to Sun (AU) (R8)
c     ierr  : Error code (I4)
c                0 : nominal run
c               10 : unknown observatory code
c              100 : date of call earlier than xjdbeg
c              200 : date of call later   than xjdend
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) code
Cf2py intent(in) t
Cf2py intent(out) pos
Cf2py intent(out) vel
Cf2py intent(out) r
Cf2py intent(out) ierr

      implicit none

      integer*4
     $  code, ierr, istat

      real*8
     $  pos(3), t, vel(3), km2AU, pos_b(3), vel_b(3), r

      parameter
     $  (km2AU = 149597870.691d0)

      ierr = 0
      if (code .eq. 1) then
         ierr = 10
         return
      else if (code .eq. 2) then
      else

c Get heliocentric position of Earth.
         call newcomb (t, pos)
         pos(1) = -pos(1)
         pos(2) = -pos(2)
         pos(3) = -pos(3)
         vel(1) = 0.d0
         vel(2) = 0.d0
         vel(3) = 0.d0

c Now move to barycenter. First get barycenter position (ecliptic).
         call BaryXV (t, pos_b, vel_b, istat)
         if (ierr .ne. 0) then
            write (6, *) 'Problem while getting barycenter position'
         end if

c Convert barycenter position to Equatorial.
         call equat_ecl (-1, pos_b, pos_b, ierr)
         if (ierr .ne. 0) then
            write (6, *) 'Problem in conversion ecliptic -> equatorial'
         end if
         call equat_ecl (-1, vel_b, vel_b, ierr)
         if (ierr .ne. 0) then
            write (6, *) 'Problem in conversion ecliptic -> equatorial'
         end if
         pos(1) = pos(1) - pos_b(1)
         pos(2) = pos(2) - pos_b(2)
         pos(3) = pos(3) - pos_b(3)
         vel(1) = vel(1) - vel_b(1)
         vel(2) = vel(2) - vel_b(2)
         vel(3) = vel(3) - vel_b(3)
         if (code .eq. 500) then
         else
            ierr = 10
            return
         end if

c Finally, computes distance from observatory to Sun.
         r = dsqrt((pos(1) + pos_b(1))**2 + (pos(2) + pos_b(2))**2
     $     + (pos(3) + pos_b(3))**2)
      end if

      end

      subroutine PlanetElem (ind, jday, a, e, inc, node, peri, capm,
     $  ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives the osculating elements in heliocentric reference
c frame of a planet at a given time. From given elements and rates.
c Valid roughly from 1900 to 2100.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     ind   : Planet index (I4)
c                1 : Mercury
c                2 : Venus
c                3 : Earth
c                4 : Mars
c                5 : Jupiter
c                6 : Saturn
c                7 : Uranus
c                8 : Neptune
c     jday  : Time of elements (Julian day) (R8)
c
c OUTPUT
c     a     : Semi-major axis (R8)
c     e     : Eccentricity (R8)
c     inc   : Inclination (R8)
c     node  : Longitude of node (R8)
c     peri  : Argument of perihelion (R8)
c     capm  : Mean anomaly (R8)
c     ierr  : Error code (I4)
c                0 : nominal run
c               10 : Unknown planet index
c               20 : jday out of range
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ind
Cf2py intent(in) jday
Cf2py intent(out) a
Cf2py intent(out) e
Cf2py intent(out) inc
Cf2py intent(out) node
Cf2py intent(out) peri
Cf2py intent(out) capm
Cf2py intent(out) ierr

      implicit none

      real*8
     $  jday, a, e, inc, node, peri, capm, Pi, TwoPi, jday_min,
     $  jday_max, drad

      integer*4
     $  ind, ierr, n_planets

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  jday_min = 2415020.0, jday_max = 2488070.0, n_planets = 8)

      real*8
     $  a_p(n_planets), e_p(n_planets), i_p(n_planets),
     $  node_p(n_planets), peri_p(n_planets), capm_p(n_planets),
     $  da_p(n_planets), de_p(n_planets), di_p(n_planets),
     $  dnode_p(n_planets), dperi_p(n_planets), dcapm_p(n_planets),
     $  jday_p, dt

c Elements are given in heliocentric reference frame.
      data
     $  jday_p /2451545.0/,
     $  a_p / 0.38709893, 0.72333199, 1.00000011, 1.52366231,
     $  5.20336301, 9.53707032, 19.19126393, 30.06896348/,
     $  e_p /0.20563069, 0.00677323, 0.01671022, 0.09341233,
     $  0.04839266, 0.05415060, 0.04716771, 0.00858587/,
     $  i_p /7.00487, 3.39471, 0.00005, 1.85061, 1.30530, 2.48446,
     $  0.76986, 1.76917/,
     $  node_p /48.33167, 76.68069, -11.26064, 49.57854, 100.55615,
     $  113.71504, 74.22988, 131.72169/,
     $  peri_p /77.45645, 131.53298, 102.94719, 336.04084, 14.75385,
     $  92.43194, 170.96424, 44.97135/,
     $  capm_p /252.25084000, 181.97973000, 100.46435000, 355.45332000,
     $  34.40438000, 49.94432000, 313.23218000, 304.88003000/,
     $  da_p/0.00000066, 0.00000092, -0.00000005, -0.00007221,
     $  0.00060737, -0.00301530, 0.00152025, -0.00125196/,
     $  de_p/0.00002527, -0.00004938, -0.00003804, 0.00011902,
     $  -0.00012880, -0.00036762, -0.00019150, 0.00002510/,
     $  di_p/-23.51, -2.86, -46.94, -25.47, -4.15, 6.11, -2.09, -3.64/,
     $  dnode_p /-446.30, -996.89, -18228.25, -1020.19, 1217.17,
     $  -1591.05, -1681.40, -151.25/,
     $  dperi_p /573.57, -108.80, 1198.28, 1560.78, 839.93, -1948.89,
     $  1312.56, -844.43/,
     $  dcapm_p /538101628.29, 210664136.06, 129597740.63, 68905103.78,
     $  10925078.35, 4401052.95, 1542547.79, 786449.21/

      ierr = 0

c Wrong planet index.
      if ((ind .lt. 1) .or. (ind .gt. n_planets)) then
         ierr = 10
         return
      end if

c Jday out of range.
      if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
         ierr = 20
         return
      end if

c Ok, do the math.
c Rates are given in arcsecond per century, and angles in degree. Also we
C have all longitudes, when we need arguments. So tranformin arguments
C and radians.
      dt = (jday - jday_p)/365.25d0/100.d0
      a = a_p(ind) + dt*da_p(ind)
      e = e_p(ind) + dt*de_p(ind)
      inc = (i_p(ind) + dt*di_p(ind)/3600.d0)*drad
      node = (node_p(ind) + dt*dnode_p(ind)/3600.d0)*drad

c Longitude of pericenter
      peri = (peri_p(ind) + dt*dperi_p(ind)/3600.d0)*drad

c Mean longitude. Change to mean anomaly
      capm = (capm_p(ind) + dt*dcapm_p(ind)/3600.d0)*drad - peri

c Now get argument of pericenter
      peri = peri - node

      return

      end

      subroutine PlanetXV (ind, jday, pos, vel, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives the position and velocity in ecliptic heliocentric
c reference frame of a planet at a given time. Uses PlanetElem.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     ind   : Planet index (I4)
c                1 : Mercury
c                2 : Venus
c                3 : Earth
c                4 : Mars
c                5 : Jupiter
c                6 : Saturn
c                7 : Uranus
c                8 : Neptune
c     jday  : Time of elements (Julian day) (R8)
c
c OUTPUT
c     pos   : Position vector (3*R8)
c     vel   : Velocity vector (3*R8)
c     ierr  : Error code (I4)
c                0 : nominal run
c               10 : Unknown planet index
c               20 : jday out of range
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ind
Cf2py intent(in) jday
Cf2py intent(out) pos
Cf2py intent(out) vel
Cf2py intent(out) ierr

      implicit none

      real*8
     $  jday, pos(3), vel(3), Pi, TwoPi, jday_min,
     $  jday_max, drad, mu

      integer*4
     $  ind, ierr, n_planets, istat

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  jday_min = 2415020.0, jday_max = 2488070.0, n_planets = 8,
     $  mu = TwoPi**2)

      real*8
     $  a, e, inc, node, peri, capm, masses(n_planets), gm

      data masses /6023600.0d0,    ! Mercury
     $              408523.7d0,    ! Venus
     $              328900.56d0,   ! Earth + Moon
     $             3098708.0d0,    ! Mars
     $                1047.349d0,  ! Jupiter
     $                3497.898d0,  ! Saturn
     $               22902.98d0,   ! Uranus
     $               19412.24d0/   ! Neptune

      ierr = 0

c Wrong planet index.
      if ((ind .lt. 1) .or. (ind .gt. n_planets)) then
         ierr = 10
         return
      end if

c Jday out of range.
      if ((jday .lt. jday_min) .or. (jday .gt. jday_max)) then
         ierr = 20
         return
      end if

c Ok, do the math.
      call PlanetElem (ind, jday, a, e, inc, node, peri, capm, istat)
      if (istat .ne. 0) then
         ierr = istat
         return
      end if
      gm = mu*(1.d0 + 1.d0/masses(ind))
      call coord_cart (gm, a, e, inc, node, peri, capm, pos(1), pos(2),
     $  pos(3), vel(1), vel(2), vel(3))

      return

      end

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine gives the cartesian coordinates of the Earth at a given
c time in a heliocentric, equatorial coordinate system. This is a self
c -contained subroutine. One can also use a progamme that reads the
c JPL ephemerides DE405 (see ss_state.f program).
c
c Tested at different times, clearly gives equatorial coordinates.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      subroutine newcomb (DJ,ROUT)
Cf2py intent(in) dj
Cf2py intent(out) rout

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision
     $  Z(8,121),Z0(80),Z1(80),Z2(80),Z3(80),Z4(80),Z5(80),Z6(80),
     $  Z7(80),Z8(80),Z9(80),ZA(80),ZB(80),ZC(8),ROU(3),CC(3,3),
     $  EL(21),C(5),RIN(3),X(3),Y(3),W(3),ROUT(3)
      integer i, k

      EQUIVALENCE(Z(1,  1),Z0(1))
      EQUIVALENCE(Z(1, 11),Z1(1))
      EQUIVALENCE(Z(1, 21),Z2(1))
      EQUIVALENCE(Z(1, 31),Z3(1))
      EQUIVALENCE(Z(1, 41),Z4(1))
      EQUIVALENCE(Z(1, 51),Z5(1))
      EQUIVALENCE(Z(1, 61),Z6(1))
      EQUIVALENCE(Z(1, 71),Z7(1))
      EQUIVALENCE(Z(1, 81),Z8(1))
      EQUIVALENCE(Z(1, 91),Z9(1))
      EQUIVALENCE(Z(1,101),ZA(1))
      EQUIVALENCE(Z(1,111),ZB(1))
      EQUIVALENCE(Z(1,121),ZC(1))
      DATA TWOPI/6.283185307179586D0/
      DATA STR/206264806.2470964D0/
      DATA RTD/57.29577951308232D0/
      DATA OBL0,C1,C2,C3/23.4522944D0,.1301250D-1,.163889D-5,.50278D-6/
      DATA TFIN/2433282.5D0/
      DATA Z0/   - 1.,  1., -   6.,-  11.,    26.,-  12.,     0.,    0.,
     1           - 1.,  2., -   3.,-   3., -   4.,    5.,     0.,    0.,
     2           - 1.,  3.,    15.,-   1., -   1.,-  18.,     0.,    0.,
     3           - 1.,  4.,    19.,-  13., -   2.,-   4.,     0.,    0.,
     4           - 1.,  0.,    33.,-  67., -  85.,-  39.,    24.,-  17.,
     5           - 1.,  1.,  2350.,-4223., -2059.,-1145., -   4.,    3.,
     6           - 1.,  2., -  65.,-  34.,    68.,-  14.,     6.,-  92.,
     7           - 1.,  3., -   3.,-   8.,    14.,-   8.,     1.,    7.,
     8           - 2.,  0., -   3.,    1.,     0.,    4.,     0.,    0.,
     9           - 2.,  1., -  99.,   60.,    84.,  136.,    23.,-   3./
      DATA Z1/   - 2.,  2., -4696., 2899.,  3588., 5815.,    10.,-   6.,
     1           - 2.,  3.,  1793.,-1735., - 595.,- 631.,    37.,-  56.,
     2           - 2.,  4.,    30.,-  33.,    40.,   33.,     5.,-  13.,
     3           - 3.,  2., -  13.,    1.,     0.,   21.,    13.,    5.,
     4           - 3.,  3., - 665.,   27.,    44., 1043.,     8.,    1.,
     5           - 3.,  4.,  1506.,- 396., - 381.,-1446.,   185.,- 100.,
     6           - 3.,  5.,   762.,- 683.,   126.,  148.,     6.,-   3.,
     7           - 3.,  6.,    12.,-  12.,    14.,   13., -   2.,    4.,
     8           - 4.,  3., -   3.,-   1.,     0.,    6.,     4.,    5.,
     9           - 4.,  4., - 188.,-  93., - 166.,  337.,     0.,    0./
      DATA Z2/   - 4.,  5., - 139.,-  38., -  51.,  189., -  31.,-   1.,
     1           - 4.,  6.,   146.,-  42., -  25.,-  91.,    12.,    0.,
     2           - 4.,  7.,     5.,-   4.,     3.,    5.,     0.,    0.,
     3           - 5.,  5., -  47.,-  69., - 134.,   93.,     0.,    0.,
     4           - 5.,  6., -  28.,-  25., -  39.,   43., -   8.,-   4.,
     5           - 5.,  7., - 119.,-  33., -  37.,  136., -  18.,-   6.,
     6           - 5.,  8.,   154.,-   1.,     0.,-  26.,     0.,    0.,
     7           - 6.,  5.,     0.,    0.,     0.,    0., -   2.,    6.,
     8           - 6.,  6., -   4.,-  38., -  80.,    8.,     0.,    0.,
     9           - 6.,  7., -   4.,-  13., -  24.,    7., -   2.,-   3./
      DATA Z3/   - 6.,  8., -   6.,-   7., -  10.,   10., -   2.,-   3.,
     1           - 6.,  9.,    14.,    3.,     3.,-  12.,     0.,    0.,
     2           - 7.,  7.,     8.,-  18., -  38.,-  17.,     0.,    0.,
     3           - 7.,  8.,     1.,-   6., -  12.,-   3.,     0.,    0.,
     4           - 7.,  9.,     1.,-   3., -   4.,    1.,     0.,    0.,
     5           - 7., 10.,     0.,    0., -   3.,    3.,     0.,    0.,
     6           - 8.,  8.,     9.,-   7., -  14.,-  19.,     0.,    0.,
     7           - 8.,  9.,     0.,    0., -   5.,-   4.,     0.,    0.,
     8           - 8., 12., -   8.,-  41., -  43.,    8., -   5.,-   9.,
     9           - 8., 13.,     0.,    0., -   9.,-   8.,     0.,    0./
      DATA Z4/   - 8., 14.,    21.,   24., -  25.,   22.,     0.,    0.,
     1           - 9.,  9.,     6.,-   1., -   2.,-  13.,     0.,    0.,
     2           - 9., 10.,     0.,    0., -   1.,-   4.,     0.,    0.,
     3           -10., 10.,     3.,    1.,     3.,-   7.,     0.,    0.,
     4             1.,- 2., -   5.,-   4., -   5.,    6.,     0.,    0.,
     5             1.,- 1., - 216.,- 167., -  92.,  119.,     0.,    0.,
     6             1.,  0., -   8.,-  47.,    27.,-   6.,     0.,    0.,
     7             2.,- 3.,    40.,-  10., -  13.,-  50.,     0.,    0.,
     8             2.,- 2.,  1960.,- 566., - 572.,-1973.,     0.,-   8.,
     9             2.,- 1., -1656.,- 616.,    64.,- 137.,     0.,    0./
      DATA Z5/     2.,  0., -  24.,   15., -  18.,-  25., -   8.,    2.,
     1             3.,- 4.,     1.,-   4., -   6.,    0.,     0.,    0.,
     2             3.,- 3.,    53.,- 118., - 154.,-  67.,     0.,    0.,
     3             3.,- 2.,   395.,- 153., -  77.,- 201.,     0.,    0.,
     4             3.,- 1.,     8.,    1.,     0.,    6.,     0.,    0.,
     5             4.,- 4.,    11.,   32.,    46.,-  17.,     0.,    0.,
     6             4.,- 3., - 131.,  482.,   460.,  125.,     7.,    1.,
     7             4.,- 2.,   525.,- 256.,    43.,   96.,     0.,    0.,
     8             4.,- 1.,     7.,-   5.,     6.,    8.,     0.,    0.,
     9             5.,- 5., -   7.,    1.,     0.,   12.,     0.,    0./
      DATA Z6/     5.,- 4.,    49.,   69.,    87.,-  62.,     0.,    0.,
     1             5.,- 3., -  38.,  200.,    87.,   17.,     0.,    0.,
     2             5.,- 2.,     3.,    1., -   1.,    3.,     0.,    0.,
     3             6.,- 6.,     0.,    0., -   4.,-   3.,     0.,    0.,
     4             6.,- 5., -  20.,-   2., -   3.,   30.,     0.,    0.,
     5             6.,- 4., - 104.,- 113., - 102.,   94.,     0.,    0.,
     6             6.,- 3., -  11.,  100., -  27.,-   4.,     0.,    0.,
     7             7.,- 6.,     3.,-   5., -   9.,-   5.,     0.,    0.,
     8             7.,- 5., -  49.,    3.,     4.,   60.,     0.,    0.,
     9             7.,- 4., -  78.,-  72., -  26.,   28.,     0.,    0./
      DATA Z7/     8.,- 7.,     1.,    3.,     5.,-   1.,     0.,    0.,
     1             8.,- 6.,     6.,-   8., -  12.,-   9.,     0.,    0.,
     2             8.,- 5.,    51.,-  10., -   8.,-  44.,     0.,    0.,
     3             8.,- 4., -  17.,-  12.,     5.,-   6.,     0.,    0.,
     4             9.,- 7.,     2.,    3.,     5.,-   3.,     0.,    0.,
     5             9.,- 6.,    13.,-  25., -  30.,-  16.,     0.,    0.,
     6             9.,- 5.,    60.,-  15., -   4.,-  17.,     0.,    0.,
     7            10.,- 7.,     2.,    5.,     7.,-   3.,     0.,    0.,
     8            10.,- 6., -   7.,   18.,    14.,    6.,     0.,    0.,
     9            10.,- 5.,     5.,-   2.,     0.,    0.,     0.,    0./
      DATA Z8/    11.,- 7.,     9.,   15.,    17.,-  10.,     0.,    0.,
     1            11.,- 6., -  12.,   42.,     8.,    3.,     0.,    0.,
     2            12.,- 7., -   4.,-   5., -   4.,    3.,     0.,    0.,
     3            13.,- 8., -  13.,-   1., -   1.,   15.,     0.,    0.,
     4            13.,- 7., -  30.,-  33., -   4.,    3.,     0.,    0.,
     5            15.,- 9.,    13.,-  16., -  17.,-  14.,     0.,    0.,
     6            15.,- 8.,   200.,-  30., -   1.,-   6.,     0.,    0.,
     7            17.,-10., -   2.,-   4., -   4.,    2.,     0.,    0.,
     8            17.,- 9., -  10.,   24.,     0.,    0.,     0.,    0.,
     9             1.,- 3., -   3.,-   1., -   2.,    5.,     0.,    0./
      DATA Z9/     1.,- 2., - 155.,-  52., -  78.,  193.,     7.,    0.,
     1             1.,- 1., -7208.,   59.,    56., 7067., -   1.,   17.,
     2             1.,  0., - 307.,-2582.,   227.,-  89.,    16.,    0.,
     3             1.,  1.,     8.,-  73.,    79.,    9.,     1.,   23.,
     4             2.,- 3.,    11.,   68.,   102.,-  17.,     0.,    0.,
     5             2.,- 2.,   136., 2728.,  4021.,- 203.,     0.,    0.,
     6             2.,- 1., - 537., 1518.,  1376.,  486.,    13.,  166.,
     7             2.,  0., -  22.,-  70., -   1.,-   8.,     0.,    0.,
     8             3.,- 4., -   5.,    2.,     3.,    8.,     0.,    0.,
     9             3.,- 3., - 162.,   27.,    43.,  278.,     0.,    0./
      DATA ZA/     3.,- 2.,    71.,  551.,   796.,- 104.,     6.,-   1.,
     1             3.,- 1., -  31.,  208.,   172.,   26.,     1.,   18.,
     2             4.,- 4., -   3.,-  16., -  29.,    5.,     0.,    0.,
     3             4.,- 3., -  43.,    9.,    13.,   73.,     0.,    0.,
     4             4.,- 2.,    17.,   78.,   110.,-  24.,     0.,    0.,
     5             4.,- 1., -   1.,   23.,    17.,    1.,     0.,    0.,
     6             5.,- 5.,     0.,    0., -   1.,-   3.,     0.,    0.,
     7             5.,- 4., -   1.,-   5., -  10.,    2.,     0.,    0.,
     8             5.,- 3., -   7.,    2.,     3.,   12.,     0.,    0.,
     9             5.,- 2.,     3.,    9.,    13.,-   4.,     0.,    0./
      DATA ZB/     1.,- 2., -   3.,   11.,    15.,    3.,     0.,    0.,
     1             1.,- 1., -  77.,  412.,   422.,   79.,     1.,    6.,
     2             1.,  0., -   3.,- 320.,     8.,-   1.,     0.,    0.,
     3             1.,  1.,     0.,-   8.,     8.,    0., -   1.,    6.,
     4             2.,- 3.,     0.,    0., -   3.,-   1.,     0.,    0.,
     5             2.,- 2.,    38.,- 101., - 152.,-  57.,     0.,    0.,
     6             2.,- 1.,    45.,- 103., - 103.,-  44.,    17., - 29.,
     7             2.,  0.,     2.,-  17.,     0.,    0.,     0.,    0.,
     8             3.,- 2.,     7.,-  20., -  30.,-  11.,     0.,    0.,
     9             3.,- 1.,     6.,-  16., -  16.,-   6.,     0.,    0./
      DATA ZC/     4.,- 2.,     1.,-   3., -   4.,-   1.,     0.,    0./

      CC(1,1) = +0.9999257180268403D0
      CC(1,2) = -0.0111782838886141D0
      CC(1,3) = -0.0048584357372842D0
      CC(2,1) = +0.0111782838782385D0
      CC(2,2) = +0.9999375206641666D0
      CC(2,3) = -0.0000271576311202D0
      CC(3,1) = +0.0048584357611567D0
      CC(3,2) = -0.0000271533600777D0
      CC(3,3) = +0.9999881973626738D0

      T  = (DJ - 2415020.D0)/36525.D0
      TP = (DJ - 2396758.D0)/365.25D0
      T18= (DJ - 2378496.D0)/36525.D0
      V = (DJ - 2451545.D0)/36525.D0
      V2 = V*V
      V3 = V*V2
      T50= (DJ - 2433282.423D0)/36525.D0
      EL( 2) =  0.16750401D-01 - T*(0.41879D-04 + 0.0503D-06*T)
      EL( 5) =  0.4908229466868880D+01 + (0.3000526416797348D-01
     1        +(0.7902463002085429D-05 +  0.5817764173314427D-07*T)*T)*T
     2        -(5.066D3 + 5.880D3*T + 0.476D3*T*T)/STR
      EL( 9) =  0.4881627934111871D+01 + (0.6283319509909086D+03
     1        + 0.5279620987282842D-05*T)*T
     2        +(0.33211D3 - 0.05299D3*T + 0.03959D3*T*T)/STR
      EL(10) =  0.6256583580497099D+01 + (0.6283019457267408D+03
     1        -(0.2617993877991491D-05 +  0.5817764173314427D-07*T)*T)*T
     2        +(5.358D3 + 5.747D3*T + 0.516D3*T*T)/STR
      EL(14) =  0.5859485105530519D+01 + (0.8399709200339593D+04
     1        +(0.4619304753611657D-04 +  0.6544984694978728D-07*T18)
     2        *T18)*T18
      EL(17) =  0.5807638839584459D+00 - (0.3375728238223387D+02
     1        -(0.3964806284113784D-04 +  0.3470781143063166D-07*T18)
     2        *T18)*T18
      EL(21) =  0.3933938487925745D+01 + (0.7101833330913285D+02
     1        -(0.1752456013106637D-03 -  0.1773109075921903D-06*T18)
     2        *T18)*T18
      EL(15) = EL(14) - EL(21)
      EL(16) = EL(14) - EL(17)
      D      = EL(14) - EL(9)
      ARG = D
      DL =    + 6466.D0*DSIN(ARG) + 13.D0*DSIN(3.*ARG)
      DR =    + 13384.D0*DCOS(ARG) +  30.D0*DCOS(3.*ARG)
      DBLARG = D + EL(15)
      ARG = DBLARG
      DL = DL +  177.D0*DSIN(ARG)
      DR = DR +  370.D0*DCOS(ARG)
      DBLARG = D - EL(15)
      ARG = DBLARG
      DL = DL -  425.D0*DSIN(ARG)
      DR = DR - 1332.D0*DCOS(ARG)
      DBLARG = 3.D0*D - EL(15)
      ARG = DBLARG
      DL = DL +   39.D0*DSIN(ARG)
      DR = DR +   80.D0*DCOS(ARG)
      DBLARG = D + EL(10)
      ARG = DBLARG
      DL = DL -   64.D0*DSIN(ARG)
      DR = DR -  140.D0*DCOS(ARG)
      DBLARG = D - EL(10)
      ARG = DBLARG
      DL = DL +  172.D0*DSIN(ARG)
      DR = DR + 360.D0*DCOS(ARG)
      ARG = D - EL(10) - EL(15)
      DL = DL - 13.D0*DSIN(ARG)
      DR = DR - 30.D0*DCOS(ARG)
      ARG = 2.D0*(EL(9)-EL(17))
      DL = DL - 13.D0*DSIN(ARG)
      DR = DR + 30.D0*DCOS(ARG)
      ARG = EL(16)
      DB =    + 577.D0*DSIN(ARG)
      DBLARG = EL(16) + EL(15)
      ARG = DBLARG
      DB = DB +  16.D0*DSIN(ARG)
      DBLARG = EL(16) - EL(15)
      ARG = DBLARG
      DB = DB -  47.D0*DSIN(ARG)
      DBLARG = EL(16) - 2.D0*(EL(9) - EL(17))
      ARG = DBLARG
      DB = DB +  21.D0*DSIN(ARG)
      ARG = ARG - EL(15)
      DB = DB + 5.D0*DSIN(ARG)
      ARG = EL(10) + EL(16)
      DB = DB + 5.D0*DSIN(ARG)
      ARG = EL(10) - EL(16)
      DB = DB + 5.D0*DSIN(ARG)
      GME = 0.43296383D+01 + .2608784648D+02*TP
      GV  = 0.19984020D+01 + .1021322923D+02*TP
      GMA = 0.19173489D+01 + .3340556174D+01*TP
      GJ  = 0.25836283D+01 + .5296346478D+00*TP
      GS  = 0.49692316D+01 + .2132432808D+00*TP
      GJ = GJ + 0.579904067D-02 * DSIN(5.D0*GS - 2.D0*GJ + 1.1719644977
     1     D0 - 0.397401726D-03 * TP)
      DG = 266.*DSIN(.555015D0 + 2.076942D0*T)
     1    +6400.*DSIN(4.035027D0 + .3525565D0*T)
     2    +(1882.D0-16.*T)*DSIN(.9990265D0 + 2.622706D0*T)
      EL(10) = DG/STR + EL(10)
      EL(12) =   DSIN(     EL(10)) * (6909794. - (17274. + 21.*T)*T)
     1         + DSIN(2.D0*EL(10)) * (  72334. -    362.*T)
     2         + DSIN(3.D0*EL(10)) * (   1050. -      8.*T)
     3         + DSIN(4.D0*EL(10)) *       17.
      RO     =                          30594. -    152.*T
     1         - DCOS(     EL(10)) * (7273841. - (18182. + 22.*T)*T)
     2         - DCOS(2.D0*EL(10)) * (  91374. -    457.*T)
     3         - DCOS(3.D0*EL(10)) * (   1446. -     11.*T)
     4         - DCOS(4.D0*EL(10)) *       25.
      EL(11) = 10.D0**(RO*1.D-09)
      DO 10 K=1,4
      DBLARG = Z(1,K)*GME + Z(2,K)*EL(10)
      ARG = DBLARG
      CS  = DCOS(ARG)
      SS  = DSIN(ARG)
      DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
      DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
   10 CONTINUE
      DO 20 K=5,44
      DBLARG = Z(1,K)*GV + Z(2,K)*EL(10)
      ARG = DBLARG
      CS  = DCOS(ARG)
      SS  = DSIN(ARG)
      DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
      DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
      DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
   20 CONTINUE
      DO 30 K=45,89
      DBLARG = Z(1,K)*GMA + Z(2,K)*EL(10)
      ARG = DBLARG
      CS  = DCOS(ARG)
      SS  = DSIN(ARG)
      DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
      DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
      DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
   30 CONTINUE
      DO 40 K=90,110
      DBLARG = Z(1,K)*GJ + Z(2,K)*EL(10)
      ARG = DBLARG
      CS  = DCOS(ARG)
      SS  = DSIN(ARG)
      DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
      DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
      DB  =(Z(7,K)*CS  + Z(8,K)*SS) +DB
   40 CONTINUE
      DO 50 K=111,121
      DBLARG = Z(1,K)*GS + Z(2,K)*EL(10)
      ARG = DBLARG
      CS  = DCOS(ARG)
      SS  = DSIN(ARG)
      DL  =(Z(3,K)*CS  + Z(4,K)*SS )+ DL
      DR  =(Z(5,K)*CS  + Z(6,K)*SS )+ DR
      DB  =(Z(7,K)*CS  + Z(8,K)*SS )+ DB
   50 CONTINUE
      C(1) = EL(11)*10.D0**(DR*1.D-09)
      C(4) = (DL + DG + EL(12))/STR + EL(9)
      C(5) = DB/STR
      C(2) = C(4)*RTD
      C(3) = C(5)*RTD
      BS = (-0.03903D3 + 0.03450D3*T - 0.00077D3*T*T)/STR
      BC = (+0.03363D3 + 0.01058D3*T + 0.00286D3*T*T)/STR
      BR = -0.01914D3/STR
      C(5) = C(5) + BS*DSIN(C(4)) + BC*DCOS(C(4)) +BR*DSIN(C(4)-EL(17))
      RS = (+0.01015D3 + 0.00035D3*T + 0.00023D3*T*T)/STR
      RC = (-0.02652D3 - 0.00002D3*T + 0.00054D3*T*T)/STR
      C(1) = C(1) + RS*DSIN(EL(10)) - RC*DCOS(EL(10))
      T2 = T*T
      T3 = T2*T
      OBL = OBL0 - C1*T - C2*T2 + C3*T3
      OBL = OBL/RTD
      RIN(1) = C(1)*DCOS(C(4))*DCOS(C(5))
      RIN(2) = C(1)*(DCOS(C(5))*DSIN(C(4))*DCOS(OBL)-DSIN(C(5))
     1*DSIN(OBL))
      RIN(3) = C(1)*(DCOS(C(5))*DSIN(C(4))*DSIN(OBL)+DSIN(C(5))
     1*DCOS(OBL))
      RDPSC=4.848136811D-06
      DPTTY=365242.2D0
      SDL=DJ
      IF (DJ-TFIN .lt. 0.) goto 101
      IF (DJ-TFIN .eq. 0.) goto 101
      IF (DJ-TFIN .gt. 0.) goto 100
  100 CONTINUE
      SDL=TFIN
  101 CONTINUE
      DT=(SDL-2415020.5D0)/DPTTY
      DT2=DT**2
      DY=DABS(DJ-TFIN)/DPTTY
      DY2=DY**2
      DY3=DY**3
      AK=((23042.53D0+139.73D0*DT+.06D0*DT2)*DY+(30.23D0-.27D0*DT)*DY2+1
     18.00D0*DY3)*RDPSC
      AW=((23042.53D0+139.73D0*DT+.06D0*DT2)*DY+(109.50D0+.39D0*DT)*DY2+
     118.32D0*DY3)*RDPSC
      BN=((20046.85D0-85.33D0*DT-.37D0*DT2)*DY+(-42.67D0-.37D0*DT)*DY2-4
     11.80D0*DY3)*RDPSC
      SINK=DSIN(AK)
      COSK=DCOS(AK)
      SINW=DSIN(AW)
      COSW=DCOS(AW)
      SINN=DSIN(BN)
      COSN=DCOS(BN)
      X(1) =-SINK*SINW   +COSK*COSW*COSN
      Y(1) =-COSK*SINW   -SINK*COSW*COSN
      W(1) =-COSW*SINN
      X(2) =SINK*COSW    +COSK*SINW*COSN
      Y(2) =COSK*COSW    -SINK*SINW*COSN
      W(2) =-SINW*SINN
      X(3) =COSK*SINN
      Y(3) =-SINK*SINN
      W(3) =COSN
      ROU(1) = 0.D0
      ROU(2) = 0.D0
      ROU(3) = 0.D0
      IF (DJ-TFIN .lt. 0.) goto 1
      IF (DJ-TFIN .eq. 0.) goto 1
      IF (DJ-TFIN .gt. 0.) goto 2
    1 CONTINUE
      DO  11  I=1,3
      ROU(I)=X(I)*RIN(1) + Y(I)*RIN(2)+W(I)*RIN(3)
   11 CONTINUE
      GO TO 3
    2 CONTINUE
      ROU(1)=      X(1)*RIN(1)   +X(2)*RIN(2)   +X(3)*RIN(3)
      ROU(2)=      Y(1)*RIN(1)   +Y(2)*RIN(2)   +Y(3)*RIN(3)
      ROU(3)=      W(1)*RIN(1)   +W(2)*RIN(2)   +W(3)*RIN(3)
    3 CONTINUE
      RIH = DSQRT(ROU(1)*ROU(1)+ROU(2)*ROU(2)+ROU(3)*ROU(3))
      RIG = DATAN2(ROU(2),ROU(1))
      DIG = DATAN2(ROU(3),DSQRT(ROU(1)*ROU(1)+ROU(2)*ROU(2)))
      RIG = RIG + 0.0D3/STR
      ROU(1) = RIH*DCOS(RIG)*DCOS(DIG)
      ROU(2) = RIH*DSIN(RIG)*DCOS(DIG)
      DO 4 I = 1,3
      ROUT(I) = CC(I,1)*ROU(1)+CC(I,2)*ROU(2)+CC(I,3)*ROU(3)
    4 CONTINUE
      RETURN
      END

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File Polygon-lib.f
c
c J.-M. Petit  Observatoire de Besançon
c Version 1 :  February 2016
c
c The purpose of this library is to provide polygon-oriented routines.
c The first and most important one is:
c     point_in_polygon (p, poly, n)
c which tels if the point "p" is inside, outside or touching the
c polygon "poly".
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      integer*4 function point_in_polygon(p, poly, n)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function checks if point "p" is inside the polygon "poly" using
c the quadrant method.
c
c This moves the given point to the origin and shifts the polygon
c accordingly. Then for each edge of the polygon, calc_walk_summand is
c called. If the sum of all returned values from these calls is +4 or
c -4, the point lies indeed inside the polygon. Otherwise, if a
c PolygonsTouching exception has been caught, the point lies on one of
c he edges of the polygon.
c
c Returns the number of nodes of the polygon, if the point lies inside,
c otherwise 1 if the point lies on the polygon and if not, 0.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     p     : Point (array (x,y)) (2*R8)
c     poly  : Array of points ((2,n+1)*R8)
c     n     : Number of vertices in the polygon (I4)
c
c OUTPUT
c     point_in_polygon: result of the call (I4)
c                Point inside polygon : n
c                Point on polygon     : 1
c                Point outside polygon: 0
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) p
Cf2py intent(in) poly
Cf2py intent(in) n

      implicit none

      integer*4
     $  n, n_max, i, walk_sum, walk, calc_walk_summand

      parameter
     $  (n_max = 100)

      real*8
     $  p(2), poly(2,n+1), moved(2,n_max+1)

      external
     $  calc_walk_summand

c Move point to origin
      do i = 1, n+1
         moved(1,i) = poly(1,i) - p(1)
         moved(2,i) = poly(2,i) - p(2)
      end do

c Computing the running sum
      walk_sum = 0
      do i = 1, n
         walk = calc_walk_summand(moved(1,i), moved(1,i+1))
         if (walk .eq. -100) then
c Point is touching the polygon
            point_in_polygon = 1
            return
         end if
c Point is not on polygon
         walk_sum = walk_sum + walk
      end do

c Final check
      if (abs(walk_sum) .eq. 4) then
         point_in_polygon = n
      else
         point_in_polygon = 0
      end if
      return

      end

      integer*4 function calc_walk_summand(p1, p2)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function calculates the summand along one edge depending on axis
c crossings.
c
c Follows the edge between two points and checks if one or both axes
c are being crossed. If They are crossed in clockwise sense, it returns
c +1 otherwise -1. Going through the origin raises the PolygonsTouching
c exception (returns -100).
c
c Returns one of -2, -1, 0, +1, +2 or raises PolygonsTouching
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     p1    : First point (array (x,y)) of edge (2*R8)
c     p2    : Second point of edge (2*R8)
c
c OUTPUT
c     calc_walk_summand: result of the call (I4)
c                Clockwise crossing        : +1
c                No crossing               : 0
c                Counter-clockwise crossing: -1
c                Diagonal crossing         : +/-2
c                Point on edge             : -100
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) p1
Cf2py intent(in) p2

      implicit none

c Indices for better readability
      integer*4
     $  x, y, summand

      parameter
     $  (x = 1, y = 2)

      real*8
     $  p1(2), p2(2), tx, ty, x_y0, y_x0

      summand = 0
c Here, we assume the 2 points are different!
c
c Checking for vertical line
      if (p1(x) .ne. p2(x)) then
         ty = p1(x)/(p1(x) - p2(x))
      else
         ty = p1(y)/(p1(y) - p2(y))
      end if

c Checking for horizontal line
      if (p1(y) .ne. p2(y)) then
         tx = p1(y)/(p1(y) - p2(y))
      else
         tx = ty
      end if

c Compute position of axis intersection
      x_y0 = p1(x) + tx*(p2(x) - p1(x))
      y_x0 = p1(y) + ty*(p2(y) - p1(y))

c Check if crossing x axis
      if ((tx .ge. 0.d0) .and. (tx .lt. 1.d0)) then
c Check if origin on edge
         if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
            calc_walk_summand = -100
            return
         end if
         x_y0 = x_y0*(p2(y) - p1(y))
         if (x_y0 .ne. 0.d0) summand = summand + sign(1.d0, x_y0)
      end if

c Check if crossing y axis
      if ((ty .ge. 0.d0) .and. (ty .lt. 1.d0)) then
c Check if origin on edge
         if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
            calc_walk_summand = -100
            return
         end if
         y_x0 = y_x0*(p1(x) - p2(x))
         if (y_x0 .ne. 0.d0) summand = summand + sign(1.d0, y_x0)
      end if
      calc_walk_summand = summand
      return

      end

      subroutine check_polygon(poly, n)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c NOTE: user must make sure the input array is at least (2,n+1) long.
c This subroutine will copy first point into index n+1, if last point
c not already same as first. It will also check that there are no two
c consecutive points the same.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     poly  : Array of points ((2,n+1)*R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      integer*4
     $  n, i, j

      real*8
     $  poly(2,*)

      if ((poly(1,n) .eq. poly(1,1)) .and. (poly(2,n) .eq. poly(2,1)))
     $  then
         n = n - 1
      else
         poly(1,n+1) = poly(1,1)
         poly(2,n+1) = poly(2,1)
      end if
      j = 1
 1000 continue
      if ((poly(1,j) .eq. poly(1,j+1))
     $  .and. (poly(2,j) .eq. poly(2,j+1))) then
         do i = j+1, n
            poly(1,i) = poly(1,i+1)
            poly(2,i) = poly(2,i+1)
         end do
         n = n - 1
         j = j - 1
      end if
      j = j + 1
      if (j .le. n) goto 1000
      return

      end

c \subroutine{coord\_cart}

      subroutine coord_cart (mu, a, ecc, inc, capo, smallo, capm,
     $  x, y, z, vx, vy, vz)

c This routine transforms delaunay variables into cartisian
c variables.
c
c ANGLES ARE GIVEN IN RADIAN !!!!
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|x, y, z, vx, vy, vz| = cartesian variables: $X, Y, Z, Px, Py, Pz$
c \end{verse}
c
c \subsubsection{Declarations}
c
Cf2py intent(in) mu
Cf2py intent(in) a
Cf2py intent(in) ecc
Cf2py intent(in) inc
Cf2py intent(in) capo
Cf2py intent(in) smallo
Cf2py intent(in) capm
Cf2py intent(out) x
Cf2py intent(out) y
Cf2py intent(out) z
Cf2py intent(out) vx
Cf2py intent(out) vy
Cf2py intent(out) vz

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, mu, x, y, z, vx, vy, vz

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
c      $i$ \\
c \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|mat| = rotation matrix \\
c \verb|q_vec, qp| = $q$ and $\dot q$ \\
c \verb|tmp| = temporary variable
c \end{verse}
c
c \subsubsection{Declarations}

      integer*4
     $  i

      real*8
     $  delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2), qp(2),
     $  sin_e, sin_i, tmp, signe, Pi, TwoPi, f, de, fp, fpp, fppp

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

c Computation of sinus and cosines of angles.

      signe = 1.d0
      if (a .lt. 0.d0) signe = -1.d0
      cos_i = dcos(inc)
      sin_i = dsqrt(1.d0 - cos_i**2)
      delau(2) = dcos(smallo)
      delau(3) = dsin(smallo)
      delau(4) = dcos(capo)
      delau(5) = dsin(capo)
      delau(1) = capm - int(capm/TwoPi)*TwoPi
      delau(6) = signe*dsqrt(mu*a*signe)
      delau(7) = abs(delau(6))*dsqrt((1.d0 - ecc**2)*signe)

c Rotation matrix.
c The rotation matrix is the composition of 3 matrices: $R_{xq} =
c  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
c \begin{displaymath}
c     R_{xq} = \left(\matrix{
c       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
c       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
c       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
c       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
c       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
c       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
c       \frac{H}{G} \cr}\right);
c \end{displaymath}

      mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
      mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
      mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
      mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
      mat(3,1) = sin_i*delau(3)
      mat(3,2) = sin_i*delau(2)

c Eccentric anomaly.
c We solve iteratively the equation:
c \begin{displaymath}
c     E - e \sin(E) = l
c \end{displaymath}
c using the accelerated Newton's method (see Danby).

      e = delau(1) + sign(.85d0, dsin(delau(1)))*ecc
      i = 0
 1000 continue
      sin_e = ecc*dsin(e)
      f = e - sin_e - delau(1)
      if (dabs(f) .gt. 1.d-14) then
         cos_e = ecc*dcos(e)
         fp = 1.d0 - cos_e
         fpp = sin_e
         fppp = cos_e
         de = -f/fp
         de = -f/(fp + de*fpp/2.d0)
         de = -f/(fp + de*fpp/2.d0 + de*de*fppp/6.d0)
         e = e + de
         i = i + 1
         if (i .lt. 20) goto 1000
         write (6, *) 'COORD_CART: No convergence: ', i, f
         write (6, *) mu, e, de
         write (6, *) a, ecc, inc
         write (6, *) capo, smallo, capm
         stop
      end if
 1100 continue

c Coordinates relative to the orbit.
c The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
c and $\vec P = R_{xq} \dot{\vec q}$, where:
c \begin{eqnarray*}
c     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
c                        \frac{GL}{\mu}\sin(E), 0\right), \\
c     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
c                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
c \end{eqnarray*}

      cos_e = dcos(e)
      sin_e = dsin(e)
      q_vec(1) = delau(6)**2*(cos_e - ecc)/mu
      q_vec(2) = delau(7)*delau(6)*sin_e/mu
      tmp = mu/(delau(6)*(1.d0 - ecc*cos_e))
      qp(1) = -sin_e*tmp
      qp(2) = delau(7)*cos_e*tmp/delau(6)

c Cartisian coordinates

      x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
      y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
      z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)
      vx = mat(1,1)*qp(1) + mat(1,2)*qp(2)
      vy = mat(2,1)*qp(1) + mat(2,2)*qp(2)
      vz = mat(3,1)*qp(1) + mat(3,2)*qp(2)

      return
      end

c \subroutine{pos\_cart}

      subroutine pos_cart (a, ecc, inc, capo, smallo, capm,
     $  x, y, z)

c This routine transforms delaunay variables into cartisian
c variables, positions only.
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|x, y, z| = cartesian variables: $X, Y, Z$
c \end{verse}
c
c \subsubsection{Declarations}
c
Cf2py intent(in) a
Cf2py intent(in) ecc
Cf2py intent(in) inc
Cf2py intent(in) capo
Cf2py intent(in) smallo
Cf2py intent(in) capm
Cf2py intent(out) x
Cf2py intent(out) y
Cf2py intent(out) z

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, x, y, z

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
c      $i$ \\
c \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|mat| = rotation matrix \\
c \verb|q_vec| = $q$ \\
c \end{verse}
c
c \subsubsection{Declarations}

      integer*4
     $  i

      real*8
     $  delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2),
     $  sin_e, sin_i, signe, Pi, TwoPi, f, de, fp, fpp, fppp

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

c Computation of sinus and cosines of angles.

      signe = 1.d0
      if (a .lt. 0.d0) signe = -1.d0
      cos_i = dcos(inc)
      sin_i = dsqrt(1.d0 - cos_i**2)
      delau(2) = dcos(smallo)
      delau(3) = dsin(smallo)
      delau(4) = dcos(capo)
      delau(5) = dsin(capo)
      delau(1) = capm - int(capm/TwoPi)*TwoPi
      delau(6) = signe*dsqrt(a*signe)
      delau(7) = abs(delau(6))*dsqrt((1.d0 - ecc**2)*signe)

c Rotation matrix.
c The rotation matrix is the composition of 3 matrices: $R_{xq} =
c  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
c \begin{displaymath}
c     R_{xq} = \left(\matrix{
c       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
c       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
c       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
c       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
c       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
c       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
c       \frac{H}{G} \cr}\right);
c \end{displaymath}

      mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
      mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
      mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
      mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
      mat(3,1) = sin_i*delau(3)
      mat(3,2) = sin_i*delau(2)

c Eccentric anomaly.
c We solve iteratively the equation:
c \begin{displaymath}
c     E - e \sin(E) = l
c \end{displaymath}
c using the accelerated Newton's method (see Danby).

      e = delau(1) + sign(.85d0, dsin(delau(1)))*ecc
      i = 0
 1000 continue
      sin_e = ecc*dsin(e)
      f = e - sin_e - delau(1)
      if (dabs(f) .gt. 1.d-14) then
         cos_e = ecc*dcos(e)
         fp = 1.d0 - cos_e
         fpp = sin_e
         fppp = cos_e
         de = -f/fp
         de = -f/(fp + de*fpp/2.d0)
         de = -f/(fp + de*fpp/2.d0 + de*de*fppp/6.d0)
         e = e + de
         i = i + 1
         if (i .lt. 20) goto 1000
         write (6, *) 'POS_CART: No convergence: ', i, f
         write (6, *) e, de
         write (6, *) a, ecc, inc
         write (6, *) capo, smallo, capm
         stop
      end if
 1100 continue

c Coordinates relative to the orbit.
c The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
c and $\vec P = R_{xq} \dot{\vec q}$, where:
c \begin{eqnarray*}
c     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
c                        \frac{GL}{\mu}\sin(E), 0\right), \\
c     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
c                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
c \end{eqnarray*}

      cos_e = dcos(e)
      sin_e = dsin(e)
      q_vec(1) = delau(6)**2*(cos_e - ecc)
      q_vec(2) = delau(7)*delau(6)*sin_e

c Cartisian coordinates

      x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
      y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
      z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)

      return
      end

c \subroutine{coord\_cart}

      subroutine PQ_cart (inc, capo, smallo, P, Q, R)

c This routine transforms delaunay variables into cartisian
c variables.
c
c ANGLES ARE GIVEN IN RADIAN !!!!
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|inc, capo, smallo| = osculating elements \\
c \verb|P, Q, R| = $\vec P$, $\vec Q$ and $\vec P \times \vec Q$ vectors
c \end{verse}
c
c \subsubsection{Declarations}
c
Cf2py intent(in) inc
Cf2py intent(in) capo
Cf2py intent(in) smallo
Cf2py intent(out) P
Cf2py intent(out) Q
Cf2py intent(out) R

      implicit none

      real*8
     $  inc, capo, smallo, P(3), Q(3), R(3)

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_i, sin_i| = sinus and cosines of $i$ \\
c \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$
c \end{verse}
c
c \subsubsection{Declarations}

      integer*4
     $  i

      real*8
     $  delau(8), cos_i, sin_i

c Computation of sinus and cosines of angles.

      cos_i = dcos(inc)
      sin_i = dsqrt(1.d0 - cos_i**2)
      delau(2) = dcos(smallo)
      delau(3) = dsin(smallo)
      delau(4) = dcos(capo)
      delau(5) = dsin(capo)

c Rotation matrix.
c The rotation matrix is the composition of 3 matrices: $R_{xq} =
c  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
c \begin{displaymath}
c     R_{xq} = \left(\matrix{
c       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
c       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
c       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
c       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
c       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
c       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
c       \frac{H}{G} \cr}\right);
c \end{displaymath}

      P(1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
      Q(1) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
      R(1) = sin_i*delau(5)
      P(2) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
      Q(2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
      R(2) = -sin_i*delau(4)
      P(3) = sin_i*delau(3)
      Q(3) = sin_i*delau(2)
      R(3) = cos_i

      return
      end

c \subroutine{osc\_el}

      subroutine osc_el (mu, x, y, z, vx, vy, vz,
     $  a, ecc, inc, capo, smallo, capm)

c This routine transforms cartisian variables into delaunay
c variables.
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|cart| = cartesian variables: $X, Y, Z, Px, Py, Pz$
c \end{verse}
c
c \subsubsection{Declarations}
Cf2py intent(in) mu
Cf2py intent(in) x
Cf2py intent(in) y
Cf2py intent(in) z
Cf2py intent(in) vx
Cf2py intent(in) vy
Cf2py intent(in) vz
Cf2py intent(out) a
Cf2py intent(out) ecc
Cf2py intent(out) inc
Cf2py intent(out) capo
Cf2py intent(out) smallo
Cf2py intent(out) capm

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, mu, x, y, z, vx, vy, vz

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_i, sin_i| = sinus and cosines of $i$ \\
c \verb|delau| = Delaunay variables:  $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|f| = $f$ true anomaly \\
c \verb|g| = $g$ argument of pericenter \\
c \verb|h_vec| = $\vec h = \vec X \times \vec P$ \\
c \verb|p_vec| = $\vec p = -\mu \frac{\vec X}{r}
c      - \vec h \times \vec P$ \\
c \verb|r| = radial distance \\
c \verb|tmp1, tmp2| = temporary variables \\
c \verb|v2| = velocity squared
c \end{verse}
c
c \subsubsection{Declarations}

      real*8
     $  delau(8), e, f, h_vec(3), p_vec(3),
     $  r, tmp1, tmp2, v2, signe, cart(6)

c Computation of angular momentum and eccentricity vector.
c \begin{eqnarray*}
c     \vec h & = & \vec X \times \vec P, \\
c     \vec p & = & -\mu \frac{\vec X}{|\vec X|}
c                  - \vec h \times \vec P
c \end{eqnarray*}

      cart(1) = x
      cart(2) = y
      cart(3) = z
      cart(4) = vx
      cart(5) = vy
      cart(6) = vz
      h_vec(1) = cart(2)*cart(6) - cart(3)*cart(5)
      h_vec(2) = cart(3)*cart(4) - cart(1)*cart(6)
      h_vec(3) = cart(1)*cart(5) - cart(2)*cart(4)
      r = 1.d0/dsqrt(cart(1)**2 + cart(2)**2 + cart(3)**2)
      p_vec(1) = -mu*cart(1)*r - h_vec(2)*cart(6) + h_vec(3)*cart(5)
      p_vec(2) = -mu*cart(2)*r - h_vec(3)*cart(4) + h_vec(1)*cart(6)
      p_vec(3) = -mu*cart(3)*r - h_vec(1)*cart(5) + h_vec(2)
     $  *cart(4)

c Computation of momenta.
c \begin{eqnarray*}
c     L & = & \mu\sqrt{\frac{1}
c                      {\frac{2\mu}{|\vec X|}-|\vec P|^2}}, \\
c     G & = & |\vec h|, \\
c     H & = &  h_z
c \end{eqnarray*}

      v2 = cart(4)**2 + cart(5)**2 + cart(6)**2
      tmp1 = 2.d0*mu*r - v2
      signe = 1.d0
      if (tmp1 .lt. 0.d0) signe = -1.d0
      delau(6) = signe*mu/dsqrt(signe*tmp1)
      delau(7) = dsqrt(h_vec(1)**2 + h_vec(2)**2 + h_vec(3)**2)
      delau(8) = h_vec(3)

      if ((cart(3) .eq. 0.d0) .and. (cart(6) .eq. 0.d0)) then
         delau(4) = 1.d0
         delau(5) = 0.d0
         tmp1 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2)
         delau(2) = p_vec(1)*tmp1
         delau(3) = p_vec(2)*tmp1
      else

c Longitude of node.
c \begin{eqnarray*}
c     \cos(h) & = & -\frac{h_y}{\sqrt{h_x^2 + h_y^2}}, \\
c     \sin(h) & = & \frac{h_x}{\sqrt{h_x^2 + h_y^2}}
c \end{eqnarray*}

         tmp1 = 1.d0/dsqrt(h_vec(1)**2 + h_vec(2)**2)
         delau(4) = -h_vec(2)*tmp1
         delau(5) = h_vec(1)*tmp1

c Argument of pericenter.
c Let us call $\vec N$ the vector derived from $\vec h$, pointing
c to the ascending node:
c \begin{displaymath}
c     \vec N = \left(-h_y, h_x, 0\right)
c \end{displaymath}
c From this, we get:
c \begin{eqnarray*}
c     \cos(g) & = & \frac{\vec N \cdot \vec p}{|\vec N||\vec p|}, \\
c     \sin(g) & = & \frac{\vec N \times \vec p}{|\vec N||\vec p|}
c                   \cdot \frac{\vec h}{|\vec h|}
c \end{eqnarray*}

         tmp2 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2 + p_vec(3)**2
     $     )
         delau(2) = (h_vec(1)*p_vec(2) - h_vec(2)*p_vec(1))*tmp1
     $     *tmp2
         delau(3) = ((h_vec(1)**2 + h_vec(2)**2)*p_vec(3)
     $     - h_vec(3)*(h_vec(1)*p_vec(1) + h_vec(2)*p_vec(2)))
     $     *tmp1*tmp2/delau(7)
      end if

c Mean anomaly
c We define $\vec X_{orb} = R_1(i) \cdot R_3(h) \vec X$. It turns
c out that $\vec X_{orb} = (r \cos(g+f), r \sin(g+f), 0)$. Hence:
c \begin{eqnarray*}
c     \cos(g+f) & = & \cos(h) X + \sin(h) Y, \\
c     \sin(g+f) & = & \cos(i) \left(\cos(h) Y - \sin(h) X\right)
c                   + \sin(i) Z
c \end{eqnarray*}
c Furthermore, we have the relation:
c \begin{displaymath}
c     \tan(\frac{E}{2}) = \sqrt{\frac{1 - e}{1 + e}}
c                         \tan(\frac{f}{2})
c \end{displaymath}
c and finally:
c \begin{displaymath}
c     l = E - e \sin(E)
c \end{displaymath}

      ecc = dsqrt(dmax1(1.d0 - signe*(delau(7)/delau(6))**2, 0.d0))
      tmp1 = (cart(1)*p_vec(1) + cart(2)*p_vec(2) + cart(3)*p_vec(3))
      tmp2 = ((cart(3)*p_vec(2) - cart(2)*p_vec(3))*h_vec(1)
     $  + (cart(1)*p_vec(3) - cart(3)*p_vec(1))*h_vec(2)
     $  + (cart(2)*p_vec(1) - cart(1)*p_vec(2))*h_vec(3))
     $  /delau(7)
      f = datan2(tmp2, tmp1)
      e = 2.d0*datan(dsqrt(signe*(1.d0 - ecc)/(1.d0 + ecc))
     $  *dtan(f/2.d0))
      capm = e - ecc*dsin(e)
      capo = datan2(delau(5), delau(4))
      smallo = datan2(delau(3), delau(2))
      a = signe*delau(6)**2/mu
      inc = dacos(dmax1(dmin1(delau(8)/delau(7),1.d0),-1.d0))

      return
      end
      subroutine equ_ecl(epsilon, poseq, posecl)
c******************************************************
c
c    Transformation of a vector from the equatorial frame to the ecliptic frame
c
c    ANGLE IN RADIAN !!!
c
c  Author  : F. Mignard OCA/CERGA
c
c
c*** INPUT
c       epsilon : obliquity
c       poseq   : input vector in equatorial frame
c
c*** OUTPUT
c       posecl  : output vector in ecliptic frame
c
c******************************************************
      implicit none

      real*8
     $  Pi, degrad, phi, poseq(3), posecl(3), epsilon, phir

c Chapront et al. 2002 gamma to O_icrs in arcsec
      parameter
     $  (Pi = 3.141592653589793238d0, phi = 0.05542d0,
     $  degrad = Pi/180.d0)

      phir    = phi/3600d0*degrad
 
      call RotZ(-phir, poseq, posecl)

      call RotX(epsilon, posecl, posecl)

      return
      end

      subroutine ecl_equ(epsilon, posecl, poseq)
c******************************************************
c
c    Transformation of a vector from the ecliptic frame to the equatorial frame
c
c    ANGLE IN RADIAN !!!
c
c  Author  : F. Mignard OCA/CERGA
c
c
c*** INPUT
c       epsilon : obliquity
c       poseq   : input vector in ecliptic frame
c
c*** OUTPUT
c       posecl  : output vector in equatorial frame
c
c******************************************************
      implicit none

      real*8
     $  Pi, degrad, phi, poseq(3), posecl(3), epsilon, phir

c Chapront et al. 2002 gamma to O_icrs in arcsec
      parameter
     $  (Pi = 3.141592653589793238d0, phi = 0.05542d0,
     $  degrad = Pi/180.d0)

      phir    = phi/3600d0*degrad

      call RotX(-epsilon, posecl, poseq)

      call RotZ(phir, poseq, poseq)

      return
      end

      subroutine RotX (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the X axis.
c This is the same as giving the coordinates of a vector that have been
c rotating by angle -alpha around the axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(1) = pos(1)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(2) = pos(2)*ca + pos(3)*sa
      posout(3) = pos(3)*ca - pos(2)*sa

      return
      end

      subroutine RotY (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the Y axis.
c This is the same as giving the coordinates of a vector that have been
c rotating by angle -alpha around the axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(2) = pos(2)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(3) = pos(3)*ca + pos(1)*sa
      posout(1) = pos(1)*ca - pos(3)*sa

      return
      end

      subroutine RotZ (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the Z axis.
c This is the same as giving the coordinates of a vector that have been
c rotating by angle -alpha around the axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(3) = pos(3)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(1) = pos(1)*ca + pos(2)*sa
      posout(2) = pos(2)*ca - pos(1)*sa

      return
      end

      subroutine equat_ecl(ieqec,v_in,v_out,ierr)
!
!     Transformation of a vector v_in(3) from :
!        equator to ecliptic  : irot = +1
!        ecliptic to equator  : irot = -1
!        at J2000.
!     The equator is assumed to be the ICRF frame
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call equat_ecl(ieqec,vv,vv,ierr)
!
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from equator to ecliptic ; -1 from ecliptic to equator.
!     v_in     :  input vector v_in(3) in the initial frame
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
Cf2py intent(in) ieqec
Cf2py intent(in) v_in
Cf2py intent(out) v_out
Cf2py intent(out) ierr

      implicit none

      real*8
     $  epsilon, v_in(3), v_out(3), ww(3), coseps, sineps, Pi, secrad

c obliquity at J20000 in arcsec

      parameter
     $  (epsilon = 84381.41d0, Pi = 3.141592653589793238d0,
     $  secrad = Pi/180.d0/3600.d0)

      integer*4
     $  ieqec, ierr

      coseps = dcos(epsilon*secrad)
      sineps = dsin(epsilon*secrad)

c regular exit
      ierr   = 0

c to allow a call like ::  call equat_ecl(ieqec,vv,vv,ierr)
      ww(1)  = v_in(1)
      ww(2)  = v_in(2)
      ww(3)  = v_in(3)

      if (ieqec .eq. 1) then
        v_out(1) =   ww(1)
        v_out(2) =   coseps*ww(2) + sineps*ww(3)
        v_out(3) = - sineps*ww(2) + coseps*ww(3)
      else if (ieqec .eq. -1) then
        v_out(1) =   ww(1)
        v_out(2) =   coseps*ww(2) - sineps*ww(3)
        v_out(3) =   sineps*ww(2) + coseps*ww(3)
      else
c anomalous exit ieqec not allowed
        ierr     =   100
      end if

      return
      end

      subroutine invar_ecl(ieqec,v_in,v_out,ierr)
!
!     Transformation of a vector v_in(3) from :
!        invariable plane to ecliptic  : irot = +1
!        ecliptic to invariable plane  : irot = -1
!        at J2000.
!     The invariable plane is given with respect to the ecliptic plane
!     (J2000) by Burkhardt, AA, 1982:
!     inclination of invariable plane:   1° 35' 13.86" =   5713.86"
!     direction of ascending node:     107° 36' 30.8"  = 387390.8"
!
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call invar_ecl(ieqec,vv,vv,ierr)
!
!     Adapted from:
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!     J-M. Petit  UBC/CNRS
!     Version 1 : September 2007
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from invariable to ecliptic ; -1 from ecliptic to invariable.
!     v_in     :  input vector v_in(3) in the initial frame
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
Cf2py intent(in) ieqec
Cf2py intent(in) v_in
Cf2py intent(out) v_out
Cf2py intent(out) ierr

      implicit none

      real*8
     $  epsilon, v_in(3), v_out(3), ww(3), coseps, sineps, Pi, secrad,
     $  omega, cosom, sinom

c obliquity at J20000 in arcsec

      parameter
     $  (epsilon = 5713.86d0, omega = 387390.8d0,
     $  Pi = 3.141592653589793238d0, secrad = Pi/180.d0/3600.d0)

      integer*4
     $  ieqec, ierr

      coseps = dcos(epsilon*secrad)
      sineps = dsin(epsilon*secrad)
      cosom = dcos(omega*secrad)
      sinom = dsin(omega*secrad)

c regular exit
      ierr   = 0

c to allow a call like ::  call invar_ecl(ieqec,vv,vv,ierr)
      ww(1)  = v_in(1)
      ww(2)  = v_in(2)
      ww(3)  = v_in(3)

      if (ieqec .eq. 1) then
        v_out(1) =   cosom*ww(1) - sinom*(coseps*ww(2) - sineps*ww(3))
        v_out(2) =   sinom*ww(1) + cosom*(coseps*ww(2) - sineps*ww(3))
        v_out(3) =   sineps*ww(2) + coseps*ww(3)
      else if (ieqec .eq. -1) then
        v_out(1) =   cosom*ww(1) + sinom*ww(2)
        v_out(2) =   coseps*(-sinom*ww(1) + cosom*ww(2)) + sineps*ww(3)
        v_out(3) = - sineps*(-sinom*ww(1) + cosom*ww(2)) + coseps*ww(3)
      else
c anomalous exit ieqec not allowed
        ierr     =   100
      end if

      return
      end

      subroutine invar_ecl_osc(ieqec, ai, ei, ii, noi, pei, mi,
     $  ao, eo, io, noo, peo, mo, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine convert osculating elements back and forth between
c invariable plane and ecliptic plane.
c Uses invar_ecl to do the work.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ieqec
Cf2py intent(in) ai
Cf2py intent(in) ei
Cf2py intent(in) ii
Cf2py intent(in) noi
Cf2py intent(in) pei
Cf2py intent(in) mi
Cf2py intent(out) ao
Cf2py intent(out) eo
Cf2py intent(out) io
Cf2py intent(out) noo
Cf2py intent(out) peo
Cf2py intent(out) mo
Cf2py intent(out) ierr
 
      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  posi(3), poso(3), veli(3), velo(3), mu, Pi, TwoPi, drad,
     $  aid, eid, iid, noid, peid, mid

      integer*4
     $  ieqec, ierr

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  mu = TwoPi**2)

      aid = ai
      eid = ei
      iid = ii
      noid = noi
      peid = pei
      mid = mi
      call coord_cart (mu, aid, eid, iid, noid, peid, mid, posi(1),
     $  posi(2), posi(3), veli(1), veli(2), veli(3))
      call invar_ecl (ieqec, posi, poso, ierr)
      call invar_ecl (ieqec, veli, velo, ierr)
      call osc_el (mu, poso(1), poso(2), poso(3), velo(1), velo(2),
     $  velo(3), ao, eo, io, noo, peo, mo)
      call ztopi (io)
      if (io .gt. Pi) io = io - Pi
      call ztopi (noo)
      call ztopi (peo)
      call ztopi (mo)

      return
      end

      subroutine invar_ecl_inc_node(ieqec, ii, noi, io, noo, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine convert inclination and longitude of node elements back
c and forth between invariable plane and ecliptic plane.
c Uses invar_ecl_osc to do the work.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ieqec
Cf2py intent(in) ii
Cf2py intent(in) noi
Cf2py intent(out) io
Cf2py intent(out) noo
Cf2py intent(out) ierr
 
      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  iid, noid

      integer*4
     $  ieqec, ierr

      data
     $  ai /10.d0/,
     $  ei /0.2d0/,
     $  pei /0.d0/,
     $  mi /0.d0/

      iid = ii
      noid = noi
      call invar_ecl_osc (ieqec, ai, ei, iid, noid, pei, mi, ao, eo, io,
     $  noo, peo, mo, ierr)

      return
      end

      subroutine ref_ecl(ieqec,v_in, v_out, eps, om, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine converts vectors back and forth between
c given reference frame and ecliptic plane.
c Reference frame is given by Omega (om), the node longitude (rotation
c around Z axis of ecliptic) and Epsilon (eps), the inclination of the
c reference frame (rotation around node axis of ecliptic).
c
c    ANGLE IN RADIAN !!!
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!     Transformation of a vector v_in(3) from :
!        given reference frame to ecliptic  : irot = +1
!        ecliptic to given reference frame  : irot = -1
!        at J2000.
!
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call ref_ecl(ieqec,vv,vv,ierr)
!
!     Adapted from:
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!     J-M. Petit  UBC/CNRS
!     Version 1 : June 2015
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from reference frame to ecliptic ; 
!                -1 ==> from ecliptic to reference frame
!     v_in     :  input vector v_in(3) in the initial frame
!     eps      :  Epsilon, inclination of ref frame [rad]
!     om       :  Omega, longitude of node [rad]
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
Cf2py intent(in) ieqec
Cf2py intent(in) v_in
Cf2py intent(in) eps
Cf2py intent(in) om
Cf2py intent(out) v_out
Cf2py intent(out) ierr

      implicit none

      real*8
     $  eps, v_in(3), v_out(3), ww(3), coseps, sineps, Pi, secrad,
     $  om, cosom, sinom

c obliquity at J20000 in arcsec

      parameter
     $  (Pi = 3.141592653589793238d0, secrad = Pi/180.d0/3600.d0)

      integer*4
     $  ieqec, ierr

      coseps = dcos(eps)
      sineps = dsin(eps)
      cosom = dcos(om)
      sinom = dsin(om)

c regular exit
      ierr   = 0

c to allow a call like ::  call ref_ecl(ieqec,vv,vv,ierr)
      ww(1)  = v_in(1)
      ww(2)  = v_in(2)
      ww(3)  = v_in(3)

      if (ieqec .eq. 1) then
        v_out(1) =   cosom*ww(1) - sinom*(coseps*ww(2) - sineps*ww(3))
        v_out(2) =   sinom*ww(1) + cosom*(coseps*ww(2) - sineps*ww(3))
        v_out(3) =   sineps*ww(2) + coseps*ww(3)
      else if (ieqec .eq. -1) then
        v_out(1) =   cosom*ww(1) + sinom*ww(2)
        v_out(2) =   coseps*(-sinom*ww(1) + cosom*ww(2)) + sineps*ww(3)
        v_out(3) = - sineps*(-sinom*ww(1) + cosom*ww(2)) + coseps*ww(3)
      else
c anomalous exit ieqec not allowed
        ierr     =   100
      end if

      return
      end

      subroutine ref_ecl_osc(ieqec, ai, ei, ii, noi, pei, mi,
     $  ao, eo, io, noo, peo, mo, eps, om, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine converts osculating elements back and forth between
c given reference frame and ecliptic plane.
c Uses ref_ecl to do the work.
c Reference frame is given by Omega (om), the node longitude (rotation
c around Z axis of ecliptic) and Epsilon (eps), the inclination of the
c reference frame (rotation around node axis of ecliptic).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ieqec
Cf2py intent(in) ai
Cf2py intent(in) ei
Cf2py intent(in) ii
Cf2py intent(in) noi
Cf2py intent(in) pei
Cf2py intent(in) mi
Cf2py intent(in) eps
Cf2py intent(in) om
Cf2py intent(out) ao
Cf2py intent(out) eo
Cf2py intent(out) io
Cf2py intent(out) noo
Cf2py intent(out) peo
Cf2py intent(out) mo
Cf2py intent(out) ierr

      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  posi(3), poso(3), veli(3), velo(3), mu, Pi, TwoPi, drad,
     $  aid, eid, iid, noid, peid, mid, eps, om

      integer*4
     $  ieqec, ierr

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  mu = TwoPi**2)

      aid = ai
      eid = ei
      iid = ii
      noid = noi
      peid = pei
      mid = mi
      call coord_cart (mu, aid, eid, iid, noid, peid, mid, posi(1),
     $  posi(2), posi(3), veli(1), veli(2), veli(3))
      call ref_ecl (ieqec, posi, poso, eps, om, ierr)
      call ref_ecl (ieqec, veli, velo, eps, om, ierr)
      call osc_el (mu, poso(1), poso(2), poso(3), velo(1), velo(2),
     $  velo(3), ao, eo, io, noo, peo, mo)
      call ztopi (io)
      if (io .gt. Pi) io = io - Pi
      call ztopi (noo)
      call ztopi (peo)
      call ztopi (mo)

      return
      end

      subroutine ref_ecl_inc_node(ieqec, ii, noi, io, noo, eps, om,
     $  ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine converts inclination and longitude of node elements back
c and forth between given reference frame and ecliptic plane.
c Uses ref_ecl_osc to do the work.
c Reference frame is given by Omega (om), the node longitude (rotation
c around Z axis of ecliptic) and Epsilon (eps), the inclination of the
c reference frame (rotation around node axis of ecliptic).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) ieqec
Cf2py intent(in) ii
Cf2py intent(in) noi
Cf2py intent(in) eps
Cf2py intent(in) om
Cf2py intent(out) io
Cf2py intent(out) noo
Cf2py intent(out) ierr

      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  iid, noid, eps, om

      integer*4
     $  ieqec, ierr

      data
     $  ai /10.d0/,
     $  ei /0.2d0/,
     $  pei /0.d0/,
     $  mi /0.d0/

      iid = ii
      noid = noi
      call ref_ecl_osc (ieqec, ai, ei, iid, noid, pei, mi, ao, eo, io,
     $  noo, peo, mo, eps, om, ierr)

      return
      end

      subroutine ztopi (var)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This subroutine resets var to be between 0 and 2Pi.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      real*8
     $  Pi, TwoPi, var

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

 1000 continue
      if (var .gt. TwoPi) then
         var = var - TwoPi
         goto 1000
      end if
 1100 continue
      if (var .lt. 0.d0) then
         var = var + TwoPi
         goto 1100
      end if

      return
      end
