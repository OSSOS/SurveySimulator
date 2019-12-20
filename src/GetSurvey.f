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
c     sur_t : Time of pointing (n*R8) [days]
c     sur_ff: Filling factor (n*R8)
c     sur_co: Observatory code (n*I4)
c     sur_x : ICRF x coordinate of observatory (n*R8)
c     sur_y : ICRF y coordinate of observatory (n*R8)
c     sur_z : ICRF z coordinate of observatory (n*R8)
c     sur_r : Distance from observatory to Sun (n*R8)
c     sur_t2: Time of pointing 2 hours later (n*R8) [days]
c     sur_x2: ICRF x coordinate of observatory 2 hours later (n*R8)
c     sur_y2: ICRF y coordinate of observatory 2 hours later (n*R8)
c     sur_z2: ICRF z coordinate of observatory 2 hours later (n*R8)
c     sur_r2: Distance from observatory to Sun 2 hours later (n*R8)
c     sur_ef: Name of efficiency function file (n*CH80)
c     sur_nr: Number of efficiency functions per pointing (n*I4)
c     sur_rt: Rates limits for efficiency function ([rad/day]) (2,n_r_max,n*R8)
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

      include 'param.inc'

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
     $  survey*(*), sur_ef(n_sur_max)*80, eff_name*80

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
Cf2py intent(in) base_name
Cf2py intent(in) len
Cf2py intent(out) i1
Cf2py intent(out) i2
Cf2py intent(out) finished
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
c     jday  : Time of pointing (R8) [days]
c     ff    : Filling factor (R8)
c     code  : Observatory code (I4)
c     pos   : Observatory position (3*R8)
c     r     : Observatory distance to Sun (R8)
c     jday2 : Time of pointing 2 hours later (R8) [days]
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

      include 'param.inc'

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

      include 'param.inc'

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
Cf2py intent(in) ra
Cf2py intent(in) dec
Cf2py intent(out) poly
Cf2py intent(out) n_e
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
Cf2py intent(in) w
Cf2py intent(in) h
Cf2py intent(in) ra
Cf2py intent(in) dec
Cf2py intent(out) poly
Cf2py intent(out) n_e

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
Cf2py intent(in) ra
Cf2py intent(in) dec
Cf2py intent(in,out) poly
Cf2py intent(in) n_e

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
Cf2py intent(in) command
Cf2py intent(in) nwmax
Cf2py intent(out) nw
Cf2py intent(out) word
Cf2py intent(out) lw

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
Cf2py intent(in) str
Cf2py intent(out) val
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

      include 'EffUtils.f'
      include 'PosVelUtils.f'
      include 'Polygon-lib.f'
