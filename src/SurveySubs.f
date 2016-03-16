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
     $  n_sur_max, n_bin_max, n_r_max, screen, keybd, verbose,
     $  lun_s, lun_h, n_e_max

      parameter
     $  (n_sur_max = 1000, n_bin_max=100, n_r_max=10, n_e_max=41,
     $  screen = 6, keybd = 5, verbose = 9,
     $  lun_s = 13, lun_h = 6)

      real*8
     $  Pi, TwoPi, drad, eps

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0, eps = 1.d-14)

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
         call GetSurvey (surnam, lun_s, n_sur_max, n_bin_max, n_r_max,
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
                           eff_l = eta(n_bin_max, sur_rt(1,1,i_sur),
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
                              eff_l = eta(n_bin_max, sur_rt(1,1,i_sur),
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

      include 'GetSurvey.f'

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
