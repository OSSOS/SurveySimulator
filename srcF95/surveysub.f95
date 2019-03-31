module surveysub

  use datadec
  use poly_lib
  use effut
  use numutils
  use getsur
  use ioutils

contains

  subroutine Detos1 (o_m, jday, hx, color, gb, ph, period, amp, surnam, seed, &
       flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, &
       jdayp, ic, surna, h_rand)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine determines if a given object is seen by the survey
! described in the directory \verb|surnam|.
! An object is described by its ecliptic (J2000) barycentric osculating
! elements given at time \verb|jday|.
! This version uses polygons to describe the footprint of the block on
! the sky.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : January 2006
! Version 2 : May 2013
! Version 3 : March 2016
! Version 4 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file).
!             Continue looping on pointings until object is detected,
!             characterized and tracked. Don't stop at first detection.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     o_m   : orbital elements of object (t_orb_m)
!     jday  : Time of elements [JD] (R8)
!     hx    : Absolute magnitude of object in 'x' band, what ever this is (R8)
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
!     gb    : opposition surge factor G, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch jday [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     surnam: Survey directory name (CH)
!
! OUTPUT
!     seed  : Random number generator seed (I4)
!     flag  : Return flag (I4): 
!                0: not found 
!                1: found, but not tracked
!                2: found and tracked
!                3: characterized, but not tracked
!                4: characterized and tracked
!     ra    : Right ascension at detection [rad] (R8)
!     dec   : Declination at detection [rad] (R8)
!     d_ra  : Right ascension rate [rad/day] (R8)
!     d_dec : Declination rate [rad/day] (R8)
!     r     : Sun-object distance [AU] (R8)
!     delta : Earth-object distance [AU] (R8)
!     m_int : Intrinsic apparent magnitude, in x-band (R8)
!     m_rand: Averaged randomized magnitude, in detection filter (R8)
!     eff   : Actual efficiency of detection (R8)
!     isur  : Identification number of survey the object was in (I4)
!     mt    : Mean anomaly at discovery [rad] (R8)
!     jdayp : Time of discovery [JD] (R8)
!     ic    : Index of color used for survey (I4)
!     surna : Detection survey name (CH10)
!     h_rand: Absolute randomized magnitude, in detection filter (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) o_m
!f2py intent(in) jday
!f2py intent(in) hx
!f2py intent(in) color
!f2py intent(in) gb
!f2py intent(in) ph
!f2py intent(in) period
!f2py intent(in) amp
!f2py intent(in) surnam
!f2py intent(in,out) seed
!f2py intent(out) flag
!f2py intent(out) ra
!f2py intent(out) dec
!f2py intent(out) d_ra
!f2py intent(out) d_dec
!f2py intent(out) r
!f2py intent(out) delta
!f2py intent(out) m_int
!f2py intent(out) m_rand
!f2py intent(out) eff
!f2py intent(out) isur
!f2py intent(out) mt
!f2py intent(out) jdayp
!f2py intent(out) ic
!f2py intent(out) surna
!f2py intent(out) h_rand
    implicit none

    integer, parameter :: screen = 6, keybd = 5, verbose = 9, &
         lun_s = 13, lun_h = 6
    type(t_orb_m) :: o_m
    real (kind=8) :: jday, hx, color(*), gb, ph, period, amp, ra, dec, &
         d_ra, d_dec, r, delta, m_int, m_rand, eff, mt, jdayp, h_rand
    integer :: seed, flag, isur, ic
    character :: surnam*(*), surna*10

    type(t_orb_m), save :: o_ml
    type(t_obspos), save :: obspos(2)
    type(t_v3d), save :: pos, pos2
    type(t_ratecut), save :: rc
    type(t_polygon), save :: poly
    type(t_pointing), save :: points(n_sur_max)
    real (kind=8), save :: h, alpha, r2, ra2, dec2, sur_mmag(n_sur_max), &
         mag_err(6), photf(3), maglim, width, height, ra_p, dec_p, ff, &
         mag_max, mag_faint, random, track, jday_o, tmp, &
         eff_lim, r_min, r_max, ang, ang_w, track_max, track_mag, &
         track_slope, angle, rate, delta2, mag_peri, dmag, &
         p(2), ra_l, dec_l, d_ra_l, d_dec_l, r_l, delta_l, &
         m_int_l, m_rand_l, eff_l, mt_l
    integer, save :: i, filt_i, flag_l, n_sur, ierr, &
         incode, outcod, i_sur, nph
    character, save :: stra*13, stdec*13
    logical, save :: debug, newpos, rate_ok, first

    data &
         first /.true./, &
         debug /.false./, &
         eff_lim /0.4d0/

    flag = 0
    flag_l = 0

    if (first) then
       first = .false.

! Opens and reads in survey definitions
       call GetSurvey (surnam, lun_s, n_sur, points, sur_mmag, ierr)
       if (ierr .ne. 0) then
          if (ierr .eq. 100) then
             write (screen, *) &
                  'GetSurvey: reached maximum number of pointings, ', n_sur
          else if (ierr .eq. 10) then
             write (screen, *) 'Unable to open survey file in ',surnam
          else if (ierr .eq. 30) then
             goto 100
          else
             write (screen, *) 'Unknown return code in read_sur.'
          end if
          stop
       end if
100    continue
! Determine overall faintest 'x' magnitude for all surveys
       mag_faint = 0.d0
       do i_sur = 1, n_sur
! sur_mmag(i_sur) in survey's filter
! mag_max in 'x' filter
          mag_max = sur_mmag(i_sur) - color(points(i_sur)%c%f)
          if (mag_max .gt. mag_faint) mag_faint = mag_max
          if (debug) then
             write (verbose, *) i_sur, mag_max, mag_faint
          end if
       end do
       if (debug) then
          write (verbose, *) 'Faintest magnitude = ', mag_faint
       end if
!       stop
    end if

! Compute approximate maximum apparent 'x' magnitude
    r_l = o_m%a*(1.d0 - o_m%e)
    h = hx - amp*0.5d0
! mag_peri in 'x' filter
    call AppMag (r_l, r_l-1.d0, 1.d0, h, gb, alpha, mag_peri, ierr)

    if (mag_peri .le. mag_faint) then
       jday_o = -1.d30
       o_ml = o_m

! loop on surveys
       do i_sur = 1, n_sur
          obspos = points(i_sur)%o_pos
          ff = points(i_sur)%ff
          rc = points(i_sur)%c%r_cut
          track_max = points(i_sur)%c%track(1)
          track_mag = points(i_sur)%c%track(2)
          track_slope = points(i_sur)%c%track(3)
          filt_i = points(i_sur)%c%f
          mag_err = points(i_sur)%c%mag_er
          photf = points(i_sur)%c%photf
          poly = points(i_sur)%poly

! Quick and dirty trick to avoid some objects on not too faint surveys:
! drop objects that are fainter, at pericenter, and as seen from the
! Sun, than the faintest magnitude recorded for that survey, in 'x'
! band.
! mag_max in 'x' filter
          mag_max = sur_mmag(i_sur) - color(filt_i)

! Any chance this survey can see the object ?
! mag_peri in 'x' filter
          if (mag_peri .le. mag_max) then

             newpos = .false.
             if (abs(obspos(1)%jday-jday_o) .gt. 0.1d0) then
! $\Delta M = \sqrt(k^2 M_tot/a^3) \Delta t$
! where $M_tot$ is the total mass of the system, and
! $k^2 = (2 \Pi / 365.25)^2$
                o_ml%m = o_m%m + sqrt(gmb) &
                     *(twopi/(o_m%a**1.5d0*365.25d0))*(obspos(1)%jday-jday)
                o_ml%m = o_ml%m - int(o_ml%m/twopi)*twopi
                call pos_cart (o_ml, pos)
                jday_o = obspos(1)%jday
                newpos = .true.
             end if
             if (debug) then
                write (verbose, *) 'Survey: ', i_sur
                write (verbose, *) pos%x, pos%y, pos%z, jday, jday_o
                write (verbose, *) obspos(1)%pos%x, obspos(1)%pos%y, &
                     obspos(1)%pos%z, obspos(1)%jday
                write (verbose, *) obspos(2)%pos%x, obspos(2)%pos%y, &
                     obspos(2)%pos%z, obspos(2)%jday
             end if
             call DistSunEcl (obspos(1)%jday, pos, r_l)
             call RADECeclXV (pos, obspos(1)%pos, delta_l, ra_l, dec_l)
             p(1) = ra_l
             p(2) = dec_l
! Get mag in actual survey filter.
             h = hx + color(filt_i) &
                  + amp*0.5d0*sin((obspos(1)%jday-jday)/period*twopi+ph)
! mag in survey's filter
             call AppMag (r_l, delta_l, obspos(1)%r, h, gb, alpha, m_int_l, ierr)
             if (ierr .ne. 0) then
                write (screen, *) 'AppMag: something''s wrong !'
                write (screen, *) 'ierr = :', ierr
                write (screen, *) 'Survey number: ', i_sur
                stop
             end if

! Format angles for output
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
                write (verbose, '(3(f8.3, 1x), a13, 1x, a13)') &
                     o_m%m/drad, o_m%peri/drad, o_m%node/drad, stra, stdec
                write (verbose, *) ra_l/drad, dec_l/drad, m_int_l, mag_max
             end if

! Still any chance to see it (comparison in survey filter band) ?
! sur_mmag(i_sur) in survey's filter
! mag in survey's filter
             if (m_int_l .le. sur_mmag(i_sur)) then

! Is the object in the FOV ?
!
! Here we use polygons.
                ierr = point_in_polygon(p, poly)
                if (debug) then
                   write (verbose, *) 'Check for FOV.'
                   write (verbose, *) poly%n, ierr
                   do i = 1, poly%n+1
                      write (verbose, *) poly%x(i)/drad, poly%y(i)/drad
                   end do
                end if
                if (ierr .gt. 0) then

! Check for chip gaps, ..., the filling factor.
                   random = ran3(seed)
                   if (debug) then
                      write (verbose, *) &
                           'In FOV of survey. Check filling factor.'
                      write (verbose, *) random, ff
                   end if
                   if (random .le. ff) then

! Well, how is its rate of motion ? Within the rate cut or not ?
                      o_ml%m = o_m%m + (twopi/(o_m%a**1.5d0*365.25d0))*(jday_o &
                           + obspos(2)%jday - obspos(1)%jday - jday)*sqrt(gmb)
                      o_ml%m = o_ml%m - int(o_ml%m/twopi)*twopi
                      call pos_cart (o_ml, pos2)
                      call DistSunEcl (obspos(2)%jday, pos2, r2)
                      call RADECeclXV (pos2, obspos(2)%pos, delta2, ra2, dec2)
                      if (debug) then
                         write (verbose, *) 'Check for second position.'
                         write (verbose, *) o_ml%m
                         write (verbose, *) pos2%x, pos2%y, pos2%z
                         write (verbose, *) delta2, ra2/drad, dec2/drad
                      end if
                      d_ra_l = ra_l - ra2
                      if (d_ra_l .gt. Pi) d_ra_l = d_ra_l - TwoPi
                      if (d_ra_l .lt. -Pi) d_ra_l = d_ra_l + TwoPi
                      d_ra_l = d_ra_l/ &
                           (obspos(2)%jday - obspos(1)%jday)*dcos(dec_l)
                      d_dec_l = (dec2 - dec_l)/(obspos(2)%jday - obspos(1)%jday)
                      rate = dsqrt(d_ra_l**2 + d_dec_l**2)
                      angle = atan2(d_dec_l/rate, d_ra_l/rate)
                      if (angle .lt. -Pi) angle = angle + TwoPi
                      if (angle .gt. Pi) angle = angle - TwoPi
                      rate_ok = (rate .ge. rc%min) .and. (rate .le. rc%max)
                      rate_ok = rate_ok .and. &
                           (dabs(rc%angle - angle) .le. rc%hwidth)
                      if (debug) then
                         write (verbose, *) 'Check for rate.'
                         write (verbose, *) rate/drad*3600.d0/24.d0, &
                              rc%min/drad*3600.d0/24.d0, &
                              rc%max/drad*3600.d0/24.d0
                         write (verbose, *) angle/drad, rc%angle/drad, &
                              rc%hwidth/drad
                         write (verbose, *) pos2%x, pos2%y, pos2%z
                      end if
                      if (rate_ok) then

! Now check for the efficiency
                         eff_l = eta(points(i_sur)%c%eff_p, &
                              points(i_sur)%c%nr, m_int_l, rate, maglim)
                         random = ran3(seed)
                         if (debug) then
                            write (verbose, *) 'Rate OK. Check detection.'
                            write (verbose, *) random, eff_l, maglim
                         end if
                         if (random .le. eff_l) then
! Compute "measured" magnitude with 1 to 3 averaged values
                            random = ran3(seed)
                            call magran (m_int_l, mag_err, seed, tmp, dmag)
                            m_rand_l = tmp
                            if (random .gt. photf(1)) then
                               call magran (m_int_l, mag_err, seed, tmp, dmag)
                               m_rand_l = (m_rand_l + tmp)/2.d0
                            end if
                            if (random .gt. photf(1)+photf(2)) then
                               call magran (m_int_l, mag_err, seed, tmp, dmag)
                               m_rand_l = (2.d0*m_rand_l + tmp)/3.d0
                            end if
! Determine efficiency of detection for that magnitude
                            eff_l = eta(points(i_sur)%c%eff_p, &
                                 points(i_sur)%c%nr, m_rand_l, rate, maglim)
! Hurray ! We found it.
                            flag_l = 1
                            if (debug) then
                               write (verbose, *) 'Hurray ! We found it.'
                            end if

! Determine if tracked
                            random = ran3(seed)
                            track = min(track_max, &
                                 1.d0 + (m_rand_l - track_mag)*track_slope)
                            if (debug) then
                               write (verbose, *) &
                                    'Checking for track. ', random, track
                            end if
                            if (random .le. track) then
                               flag_l = 2
                            end if
! Decide if characterized or not
                            if (debug) then
                               write (verbose, *) &
                                    'Checking for characterization. ', &
                                    m_rand_l, maglim, eff_l, eff_lim
                            end if
                            if (maglim .gt. 0.d0) then
                               if (m_rand_l .le. maglim) flag_l = flag_l + 2
                            else
                               if (eff_l .ge. eff_lim) flag_l = flag_l + 2
                            end if
! Record what needs to be recorded.
                            if (flag_l .gt. flag) then
                               isur = i_sur
                               surna = points(i_sur)%efnam &
                                    (1:min(len(surna),len(points(1)%efnam)))
! Converting intrinsic magnitude to 'x' band, keeping apparent
! magnitude in discovery filter
                               ic = filt_i
                               m_int = m_int_l - color(ic)
                               m_rand = m_rand_l
                               r = r_l
                               delta = delta_l
                               call AbsMag (r, delta, obspos(1)%r, m_rand, gb, &
                                    alpha, h_rand, ierr)
                               if (ierr .ne. 0) then
                                  write (screen, *) &
                                       'AbsMag: something''s wrong !'
                                  write (screen, *) 'ierr = :', ierr
                                  write (screen, *) 'Survey number: ', i_sur
                                  stop
                               end if
                               if (debug) then
                                  write (verbose, *) 'All is good, h_rand.'
                                  write (verbose, *) r, delta
                                  write (verbose, *) m_rand, alpha, h_rand
                               end if
                               flag = flag_l
                               ra = ra_l
                               dec = dec_l
                               d_ra = d_ra_l
                               d_dec = d_dec_l
                               eff = eff_l
                               mt = o_ml%m
                               jdayp = obspos(1)%jday
                            end if
! We got it, and we know if it was tracked and/or characterized.
! Return if tracked and characterized, otherwise keep looping.
                            if (flag .ge. 4) return
                         else
!                            write (6, *) 'Low efficiency: ', a, eff_l, &
!                                 random, m_int_l, rate
                         end if
                      else
!                         write (6, *) 'Rate out of range: ', a, r_min, &
!                              r_max, rate, ang_w, dabs(ang - angle)
                      end if
                   else
!                      write (6, *) 'Falling in chip gaps: ', a, ff, random
                   end if
                end if
             else
!                write (6, *) 'Too faint for this survey: ', a, i_sur, &
!                     sur_mmag(i_sur), m_int_l, hx, filt_i, color(filt_i)
             end if
          else
!             write (6, *) 'Too faint (peri) for this survey: ', a, &
!                  i_sur, mag_max, mag_peri, hx, filt_i, color(filt_i)
          end if

! End loop on surveys
       end do
    else
!       write (6, *) 'Too faint (peri) for all surveys: ', a, &
!            mag_faint, mag_peri, hx, filt_i, color(filt_i)
    end if

    return

  end subroutine Detos1

end module surveysub
