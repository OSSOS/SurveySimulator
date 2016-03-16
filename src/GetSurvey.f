      subroutine GetSurvey (survey, lun_s, n_max, nb_max, nr_max, n_sur,
     $  sur_w, sur_h, sur_ra, sur_de, sur_t, sur_ff, sur_co, sur_x,
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
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     survey: Survey directory name (CH)
c     lun_s : Logical unit for file (I4)
c     n_max : Maximum number of pointings to read (I4)
c     nb_max: Maximum number of bin in efficiency function (I4)
c
c OUTPUT
c     n_sur : Number of pointings read (I4)
c     sur_w : Width of FOV (n*R8)
c     sur_h : Height of FOV (n*R8)
c     sur_ra: RA of pointing (n*R8)
c     sur_de: DEC of pointing (n*R8)
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
c     sur_en: Number of bins in efficiency function (nr_max,n*I4)
c     sur_eb: Bin centers for efficiency function (nb_max,nr_max,n*R8)
c     sur_em: Efficiency at bin center (nb_max,nr_max,n*R8)
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
c     sur_ml: Limiting magnitude of survey (nr_max,n*R8)
c     sur_f : Filter used for this survey (n*I4)
c     ierr  : Error code (I4)
c                0 : nominal run
c              100 : Maximum number of objects reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer*4
     $  n_bin_max, n_r_max

      parameter
     $  (n_bin_max=100, n_r_max=10)

      integer*4
     $  n_max, nb_max, nr_max

      real*8
     $  sur_w(*), sur_h(*), sur_ra(*), sur_de(*), sur_t(*),
     $  sur_ff(*), sur_x(*), sur_y(*), sur_z(*), sur_r(*),
     $  sur_t2(*), sur_x2(*), sur_y2(*), sur_z2(*), sur_r2(*),
     $  sur_eb(nb_max, nr_max, *), sur_em(nb_max, nr_max, *),
     $  width, height, ra_p, dec_p, jday_p, ff, obspos(3), ros,
     $  rates(2, n_r_max), sur_rt(2, nr_max, *),
     $  eff_b(n_bin_max, n_r_max), eff_m(n_bin_max, n_r_max),
     $  sur_mm(*), sur_rn(*), sur_rx(*), sur_an(*), sur_ph(3,*),
     $  sur_aw(*), sur_ta(*), sur_tm(*), sur_ts(*), sur_dm(6,*), 
     $  sur_ml(nr_max,*), eta, obspos2(3), ros2, maglim(n_r_max),
     $  jday_p2, rate_c(4), d_mag(6), rate, track(3), photf(3), tmp

      integer*4
     $  sur_co(*), sur_en(nr_max,*), sur_nr(*), sur_f(*), nr, j,
     $  eff_n(n_r_max), code, i, lun_s, n_sur, ierr, i1, i2, filt_i

      character
     $  survey*(*), sur_ef(*)*(*), eff_name*80

      logical
     $  finished

      external eta

c Open and read in survey definitions
      call read_file_name (survey, i1, i2, finished, len(survey))
      n_sur = 0
 200  continue
         call read_sur (survey(i1:i2), nb_max, lun_s, width, height,
     $     ra_p, dec_p, jday_p, ff, code, obspos, ros, jday_p2,
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
         sur_w(n_sur) = width
         sur_h(n_sur) = height
         sur_ra(n_sur) = ra_p
         sur_de(n_sur) = dec_p
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
                  ros = eta(nb_max, rates, sur_nr(n_sur),
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

      subroutine read_sur (dirn, nb_max, lun_in, w, h, ra, dec, jday,
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
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     dirn  : Name of directory with survey definition (CH)
c     nb_max: Maximum number of bin in efficiency function (I4)
c     lun_in: File unit (I4)
c
c OUTPUT
c     w     : field of view width (R8)
c     h     : field of view height (R8)
c     ra    : Pointing RA (R8)
c     dec   : Pointing DEC (R8)
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

      implicit none

      integer*4
     $  nw_max, nb_max

      real*8
     $  Pi, drad, TwoHours

      parameter
     $  (Pi = 3.141592653589793d0, drad = Pi/180.0D0, nw_max = 10,
     $  TwoHours = 2.d0/24.d0)

      real*8
     $  w, h, ra, dec, ff, pos(3), eff_b(nb_max, *), eff_m(nb_max, *),
     $  jday, r, vel(3), photf(3), maglim(*),
     $  pos2(3), r2, jday2, rate_c(4), track(3), d_mag(*), rates(2,*)

      integer
     $  lun_in, ierr, j, code, eff_n(*), nw, lw(nw_max), lun_e, ierr_e,
     $  i1, i2, i3, i4, filt_i, nr

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
      call read_eff (fname, nb_max, lun_e, eff_b, eff_m, eff_n,
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

      subroutine read_eff (filen, nb_max, lun_in, bin, eff, eff_n,
     $  rates, nrates, rate_c, mag_er, photf, track, maglim, filt_i,
     $  ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine opens and reads in efficiency file.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     filen : object element file name
c     nb_max: Maximum number of bin in efficiency function (I4)
c     lun_in: File unit
c
c OUTPUT
c     bin   : Magnitude bin center (nb_max,n*R8)
c     eff   : Efficiency for that bin (nb_max,n*R8)
c     eff_n : Number of bins in efficiency (n*I4)
c               -3 : piecewise linear function
c               -2 : double hyperbolic tangent
c               -1 : single hyperbolic tangent
c               >0 : number of bins in lookup table
c     rates : Rates limits for efficiency function ([rad/day]) (2,n*R8)
c     nrates: Number of efficiency fucntions (I4)
c     rate_c: rate cut parameters ([rad/day] and [rad]) (4*R8)
c     mag_er: Magnitude error parameters (6*R8)
c     photf : Photometric peasurments fractions (3*R8)
c     track : tracking fraction parameters (3*R8)
c     maglim: Limiting magnitude of survey (n*R8)
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

      implicit none

      real*8
     $  pi, drad

      integer*4
     $  nw_max, nb_max

      parameter
     $  (pi = 3.141592653589793d0, drad = pi/180.0d0, nw_max = 10)

      real*8
     $  bin(nb_max, *), eff(nb_max, *), rate_c(4), mag_er(*), track(3),
     $  rates(2,*), photf(3), maglim(*)

      integer
     $  lun_in, ierr, eff_n(*), eq_ind, nw, lw(nw_max), i, nrates, j,
     $  filt_i

      character
     $  line*100, filen*(*), word(nw_max)*80

      logical
     $  rcut, tr, fi, mag, in_rates, in_func, rate(0:10), ph

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

      include 'EffUtils.f'
      include 'PosVelUtils.f'
      include 'Polygon-lib.f'
