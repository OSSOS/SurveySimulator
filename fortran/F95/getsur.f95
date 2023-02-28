module getsur

  use datadec
  use xvutils
  use poly_lib
  use effut
  use ioutils

contains

  subroutine create_ears(ra, dec, poly)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine generates a polygon that gives the footprint of the CCD on the
! sky given the locatio of the center.
!
! Author: Jean-Marc Petit
! version 1: September 2017
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     ra    : Right ascension of center of CCD [rad] (R8)
!     dec   : Declination of center of CCD [rad] (R8)
!
! OUPUT
!     poly  : Polygons that represents the CCD
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ra
!f2py intent(in) dec
!f2py intent(out) poly
    implicit none

    type(t_polygon), intent(out) :: poly
    real (kind=8), intent(in) :: ra, dec
    real (kind=8) :: dra, h, w, e_w, e_h

! Below values are from assuming a size of 1dx1d for the "central"
! square and add 1/2 height and 1/9 width ears. However, looking at
! the resulting polygons, they seem to be somewhat too large
!      data
!     $  h /0.008726646259971648d0/,
!     $  w /0.008726646259971648d0/,
!     $  e_h /0.004363323129985824d0/,
!     $  e_w /0.0019392547244381439d0/
!
! Now, maybe we can do better. Let's look at image 1805373p.fits.
! Pointing is (13:32:29.29; -9:31:03.9) or (203.122042; -9.517750)
! Now, the corners of the "central" square are:
! (13:30:29.97; -10:00:48.4) or (202.624875; -10.013444)
! (13:34:27.66; -10:00:36.0) or (203.615250; -10.010000)
! (13:34:26.45; -09:00:58.4) or (203.610208; -09.016222)
! (13:30:29.44; -09:01:11.4) or (202.622667; -09.019833)
! This is a "square" of 0.9753d x 0.9937d
!
! Now, add the ears:
! (13:34:53.45; -09:16:26.2) or (203.722708; -09.273944)
! (13:34:27.30; -09:16:26.1) or (203.613750; -09.273917)
! (13:34:27.89; -09:45:09.2) or (203.616208; -09.752556)
! (13:34:54.05; -09:45:04.4) or (203.725208; -09.751222)
! This is a "rectangle" of 0.1075d x 0.4780d
!
! (13:30:29.09; -09:16:38.9) or (202.621208; -09.277472)
! (13:30:02.88; -09:16:42.3) or (202.512000; -09.278417)
! (13:30:03.18; -09:45:19.7) or (202.513250; -09.755472)
! (13:30:29.33; -09:45:21.5) or (202.622208; -09.755972)
! This is a "rectangle" of 0.1076d x 0.4778d
!
! So let's make the area 0.9753 x 0.9937 + 2 x 0.010755 x 0.4779
! 0.979435239 sq.deg. The CCDs are 2048x4612 with pixels
! 0.18689x0.18689" resulting in an area of 0.0254558 sq.deg. for each
! CCD, or 1.0182 sq.deg. Oops, there is a big problem. I guess it's the
! size of the pixels that's to big. Applying this size to the cenral
! square yield something reasonable, with pixels covering less than the
! whole size of the array. But the ears are in trouble. The pixel heigh
! seems to be too large. Actually, the main problem is that the WCS is
! not good enough, and the central horizontal gap has vanished, and
! some pixels are enven overlapping. I need to get an image with
! Stephen's header.
! Officially, the header says the size of the pixel is 0.185"x0.185",
! somewhat smaller than what JJ says, but I'll stick to what JJ says.
! Actually, I cannot use Stephen's header with DS9 as the latter
! doesn't know how to use the PVs, and is rather using the CDs.
!
! The horizontal size of the central square shuold be greater than
! 9x2112 = 19008 pixels or 3552.4". This implies a gap of ~32 pixels
! between the chips, in addition to the oversans, 32 pixels on each
! side. I'll assume the small horizontal gap is similar in the center
! of the frame. Then for the same 1d full size, the wisth of the large
! gaps is 327 and 328 pixels.
!
! From all this, I assume a half width of 0.5d, a half height of 0.5d,
! the width of the ears (32+2112)*0.18689" = 1/9d, and half height of
! (2*4644+32)*0.18689" = 0.483837d
!      data
!     $  h /0.008726646259971648d0/,
!     $  w /0.008726646259971648d0/,
!     $  e_h /0.004222278224995352d0/,
!     $  e_w /0.0019392547244381439d0/
!
! The above gives a lot of overlap. Using the PVs and CDs, I determined
! the exact footprint of a series of MegaPrime40 frames (in ~/Research
! /OSSOS/src/MegaPrime40.FOV) and averaged them. In addition to the
! positions below, I also determined the pixel size: 0.186" x 0.1849"
    data &
         h /0.008678d0/, &
         w /0.008545d0/, &
         e_h /0.004169d0/, &
         e_w /0.001913d0/

    poly%n = 12
    poly%y(1) = dec - h
    poly%y(12) = poly%y(1)
    poly%y(13) = poly%y(1)
    poly%y(6) = dec + h
    poly%y(7) = poly%y(6)
    poly%y(2) = dec - e_h
    poly%y(3) = poly%y(2)
    poly%y(10) = poly%y(2)
    poly%y(11) = poly%y(2)
    poly%y(4) = dec + e_h
    poly%y(5) = poly%y(4)
    poly%y(8) = poly%y(4)
    poly%y(9) = poly%y(4)

    dra = w/dcos(poly%y(1))
    poly%x(1) = ra - dra
    poly%x(13) = poly%x(1)
    poly%x(12) = ra + dra
    dra = w/dcos(poly%y(2))
    poly%x(2) = ra - dra
    poly%x(11) = ra + dra
    dra = (w + e_w)/dcos(poly%y(2))
    poly%x(3) = ra - dra
    poly%x(10) = ra + dra
    dra = (w + e_w)/dcos(poly%y(4))
    poly%x(4) = ra - dra
    poly%x(9) = ra + dra
    dra = w/dcos(poly%y(4))
    poly%x(5) = ra - dra
    poly%x(8) = ra + dra
    dra = w/dcos(poly%y(6))
    poly%x(6) = ra - dra
    poly%x(7) = ra + dra
    return
  end subroutine create_ears

  subroutine create_rectangle(w, h, ra, dec, poly)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine generates a polygon that gives the footprint of a rectangle on
! the sky given the location of the center, half-width and half-height.
!
! Author: Jean-Marc Petit
! version 1: September 2017
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     w     : Half width of rectangle [rad] (R8)
!     h     : half height of rectagle [rad] (R8)
!     ra    : Right ascension of center of rectangle [rad] (R8)
!     dec   : Declination of center of rectangle [rad] (R8)
!
! OUPUT
!     poly  : Polygons that represents the CCD
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) w
!f2py intent(in) h
!f2py intent(in) ra
!f2py intent(in) dec
!f2py intent(out) poly

    implicit none

    type(t_polygon), intent(out) :: poly
    real (kind=8), intent(in) :: w, h, ra, dec
    real (kind=8) :: dra

    poly%n = 4
    poly%y(1) = dec - h
    poly%y(4) = poly%y(1)
    poly%y(2) = dec + h
    poly%y(3) = poly%y(2)
    poly%y(5) = poly%y(1)

    dra = w/dcos(poly%y(1))
    poly%x(1) = ra - dra
    poly%x(4) = ra + dra
    dra = w/dcos(poly%y(2))
    poly%x(2) = ra - dra
    poly%x(3) = ra + dra
    poly%x(5) = poly%x(1)
    return
  end subroutine create_rectangle

  subroutine create_poly(ra, dec, poly)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine generates a polygon that gives the footprint of a polygon on
! the sky given the location of the center adn the shape of the polygon
!
! Author: Jean-Marc Petit
! version 1: September 2017
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     ra    : Right ascension of center of polygon [rad] (R8)
!     dec   : Declination of center of polygon [rad] (R8)
!
! INPUT/OUPUT
!     poly  : Polygons that represents the CCD
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ra
!f2py intent(in) dec
!f2py intent(in,out) poly

    implicit none

    type(t_polygon), intent(out) :: poly
    real (kind=8), intent(in) :: ra, dec
    real (kind=8) :: dra
    integer :: j

    do j = 1, poly%n
       poly%y(j) = dec + poly%y(j)
       poly%x(j) = ra + poly%x(j)/dcos(poly%y(j))
    end do
    poly%x(poly%n+1) = poly%x(1)
    poly%y(poly%n+1) = poly%y(1)

    return
  end subroutine create_poly

  subroutine read_eff (filen, lun_in, c, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine opens and reads in efficiency file.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : April 2013
!             Changed to read new pointings and efficiency file format
! Version 5 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     filen : object element file name
!     lun_in: File unit
!
! OUTPUT
!     c     : characterization of pointing (charact)
!     c%eff_p(j)%n : Number of bins in efficiency (I4)
!               -4 : so-called "square" function
!               -3 : piecewise linear function
!               -2 : double hyperbolic tangent
!               -1 : single hyperbolic tangent
!               >0 : number of bins in lookup table
!     c%f   : filter index (I4)
!                1 : g
!                2 : r
!                3 : i
!                4 : z
!                5 : u
!                6 : B
!                7 : V
!                8 : R
!                9 : I
!     c%r_cut : global rate cuts of search (t_ratecut)
!     c%mag_er : parameters for magnitude uncertainty and skew (6*R8)
!     c%photf : fractions of cases with 1, 2 and 3 photometric measures (3*R8)
!     c%track : tracking fraction parameters (3*R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : unable to open filen
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py inten(in) filen
!f2py inten(in) lun_in
!f2py intent(out) c
!f2py intent(out) ierr
    implicit none

    type(t_charact), intent(out) :: c
    integer, intent(in) :: lun_in
    integer, intent(out) :: ierr
    character(*), intent(in) :: filen
    integer :: eq_ind, nw, lw(nw_max), i, j
    character(256) :: line
    character(80) :: word(nw_max)
    logical rcut, tr, fi, mag, in_rates, in_func, rate(0:n_r_max), ph

    rcut = .false.
    tr = .false.
    fi = .false.
    mag = .false.
    in_rates = .false.
    in_func = .false.
    ph = .false.
    rate = .false.
    rate(0) = .true.

    ierr = 0
    open (unit=lun_in, file=filen, status='old', err=1000)
    c%nr = 0

1500 continue
    read (lun_in, '(a)', err=1500, end=3000) line
    if (line(1:1) .eq. '#') goto 1500
    eq_ind = index(line, '=')
    if (eq_ind .le. 0) goto 1500
    call parse (line(1:eq_ind-1), nw_max, nw, word, lw)
    if (nw .ne. 1) goto 1500
    if (word(1)(1:lw(1)) .eq. 'rate_cut') then
       read (line(eq_ind+1:), *, err=1500) c%r_cut%min, c%r_cut%max, &
            c%r_cut%angle, c%r_cut%hwidth
! Change rates to rad/day
       c%r_cut%min = c%r_cut%min*24.d0/3600.d0*drad
       c%r_cut%max = c%r_cut%max*24.d0/3600.d0*drad
! Change angles to radian
       c%r_cut%angle = c%r_cut%angle*drad
       c%r_cut%hwidth = c%r_cut%hwidth*drad
       rcut = .true.
       in_rates = .false.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'mag_error') then
       read (line(eq_ind+1:), *, err=1500, end=1500) (c%mag_er(i),i=1,6)
       c%mag_er(2) = log10(c%mag_er(2)/c%mag_er(1))/(c%mag_er(3)-21.d0)
       mag = .true.
       in_rates = .false.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'phot_frac') then
       read (line(eq_ind+1:), *, err=1500, end=1500) (c%photf(i),i=1,3)
       ph = .true.
       in_rates = .false.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'track_frac') then
       read (line(eq_ind+1:), *, err=1500) (c%track(i),i=1,3)
       tr = .true.
       in_rates = .false.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'filter') then
       call parse (line(eq_ind+1:), nw_max-1, nw, word(2:), lw(2:))
       if (word(2)(1:1) .eq. 'g') then
          c%f = 1
       else if (word(2)(1:1) .eq. 'r') then
          c%f = 2
       else if (word(2)(1:1) .eq. 'i') then
          c%f = 3
       else if (word(2)(1:1) .eq. 'z') then
          c%f = 4
       else if (word(2)(1:1) .eq. 'u') then
          c%f = 5
       else if (word(2)(1:1) .eq. 'B') then
          c%f = 6
       else if (word(2)(1:1) .eq. 'V') then
          c%f = 7
       else if (word(2)(1:1) .eq. 'R') then
          c%f = 8
       else if (word(2)(1:1) .eq. 'I') then
          c%f = 9
       else
          goto 1500
       end if
       fi = .true.
       in_rates = .false.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'rates') then
       j = c%nr
       if (rate(c%nr)) then
          j = j + 1
       end if
       read (line(eq_ind+1:), *, err=1500) c%eff_p(j)%min, c%eff_p(j)%max
       c%nr = j
! START comment out in production mode
!       write (18, *) 'nr = ', c%nr, c%eff_p(j)%min, c%eff_p(j)%max
! END comment out in production mode
! Change rates to rad/day
       c%eff_p(j)%min = c%eff_p(j)%min*24.d0/3600.d0*drad
       c%eff_p(j)%max = c%eff_p(j)%max*24.d0/3600.d0*drad
       c%eff_p(j)%mag_lim = -1.d0
       in_rates = .true.
       in_func = .false.
    else if (word(1)(1:lw(1)) .eq. 'function') then
       if (in_rates) then
          call parse (line(eq_ind+1:), nw_max-1, nw, word(2:), lw(2:))
          if (word(2) .eq. 'single') c%eff_p(c%nr)%n = -1
          if (word(2) .eq. 'double') c%eff_p(c%nr)%n = -2
          if (word(2) .eq. 'linear') c%eff_p(c%nr)%n = -3
          if (word(2) .eq. 'square') c%eff_p(c%nr)%n = -4
          if (word(2) .eq. 'lookup') c%eff_p(c%nr)%n = 0
          in_func = .true.
       end if
    else if (word(1)(1:lw(1)) .eq. 'linear_param') then
       if (in_func .and. (c%eff_p(c%nr)%n .eq. -3)) then
          read (line(eq_ind+1:), *, err=1500, end=1500) &
               (c%eff_p(c%nr)%e(i),i=1,3)
          rate(c%nr) = .true.
!          in_rates = .false.
          in_func = .false.
       end if
    else if (word(1)(1:lw(1)) .eq. 'single_param') then
       if (in_func .and. (c%eff_p(c%nr)%n .eq. -1)) then
          read (line(eq_ind+1:), *, err=1500, end=1500) &
               (c%eff_p(c%nr)%e(i),i=1,3)
          rate(c%nr) = .true.
!          in_rates = .false.
          in_func = .false.
       end if
    else if (word(1)(1:lw(1)) .eq. 'double_param') then
       if (in_func .and. (c%eff_p(c%nr)%n .eq. -2)) then
          read (line(eq_ind+1:), *, err=1500, end=1500) &
               (c%eff_p(c%nr)%e(i),i=1,4)
          rate(c%nr) = .true.
!          in_rates = .false.
          in_func = .false.
       end if
    else if (word(1)(1:lw(1)) .eq. 'square_param') then
       if (in_func .and. (c%eff_p(c%nr)%n .eq. -4)) then
          read (line(eq_ind+1:), *, err=1500, end=1500) &
               (c%eff_p(c%nr)%e(i),i=1,4)
          rate(c%nr) = .true.
!          in_rates = .false.
          in_func = .false.
       end if
    else if (word(1)(1:lw(1)) .eq. 'lookup_param') then
       if (in_func .and. (c%eff_p(c%nr)%n .ge. 0)) then
          i = c%eff_p(c%nr)%n + 1
          read (line(eq_ind+1:), *, err=1500, end=1500) &
               c%eff_p(c%nr)%b(i), c%eff_p(c%nr)%e(i)
          c%eff_p(c%nr)%n = i
          rate(c%nr) = .true.
       end if
    else if (word(1)(1:lw(1)) .eq. 'mag_lim') then
       if (in_rates) then
          read (line(eq_ind+1:), *, err=1500, end=1500) c%eff_p(c%nr)%mag_lim
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
       ierr = -10
       return
    end if
    if (.not. tr) then
       write (6, *) 'Survey file: ', filen
       write (6, *) 'ERROR: tracking parameters not defined.'
       ierr = -10
       return
    end if
    if (.not. fi) then
       write (6, *) 'Survey file: ', filen
       write (6, *) 'ERROR: filter not defined.'
       ierr = -10
       return
    end if
3010 continue
    if ((c%nr .gt. 0) .and. .not. rate(c%nr)) then
       c%nr = c%nr - 1
       goto 3010
    end if
    if (c%nr .le. 0) then
       write (6, *) 'Survey file: ', filen
       write (6, *) 'ERROR: no efficiency function defined.'
       ierr = -10
       return
    end if
    if (.not. mag) then
       c%mag_er(1) = 0.026d0
       c%mag_er(2) = 0.33d0
       c%mag_er(3) = 24.45d0
       c%mag_er(4) = 0.7d0
       c%mag_er(5) = 23.7d0
       c%mag_er(6) = -0.3d0
       write (6, *) 'Survey file: ', filen
       write (6, *) &
            'WARNING: magnitude error parameters not defined, using ' &
            //'default values:', c%mag_er(1), ', ', c%mag_er(2), ', ', &
            c%mag_er(3), ', ', c%mag_er(4), ', ',  c%mag_er(5), ', ', &
            c%mag_er(6)
       c%mag_er(2) = log10(c%mag_er(2)/c%mag_er(1))/(c%mag_er(3)-21.d0)
    end if
    if (.not. ph) then
       c%photf(1) = 1.d0
       c%photf(2) = 0.d0
       c%photf(3) = 0.d0
       write (6, *) 'Survey file: ', filen
       write (6, *) &
            'WARNING: photometric measurements fractions not defined, ' &
            //'using default values:', c%photf(1), ', ', c%photf(2), &
            ', ', c%photf(3)
    end if
    do j = 1, c%nr
       if (rate(j) .and. (c%eff_p(j)%mag_lim .le. 0.d0)) then
          write (6, *) 'Survey file: ', filen, '; rate range: ', &
               c%eff_p(j)%min, ' - ', c%eff_p(j)%max
          write (6, *) &
               'WARNING: limiting magnitude not defined, ' &
               //'using efficiency instead.'
       end if
    end do
    ierr = 0
    return

  end subroutine read_eff

  subroutine get_code(code_in, dirn, code_out)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! this routine checks to see if the observatory code actually references
! a file that should then contain a JPL state vector CSV file
! when a LUN is returned its assigned values starting at 501
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! JJ Kavelaars National Research Council of Canada
! Version 1 : November 2022
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     code_in  : string from pointings file that holds the observatory code
!
! OUTPUT
!     code_out : result integer code (can be observatory code or lun of open file)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) code_in
!f2py intent(out) code_out
    implicit none
    character(*), intent(IN) :: code_in, dirn
    integer, intent(OUT) :: code_out
  
    character(len=30) :: fmt, fname
    integer :: ierr, j
    integer :: vector_file_lun

    data vector_file_lun /500/
    save vector_file_lun
    
    j=len_trim(code_in)
    write(fmt, '(Ai0A)') "(I",j,")"
    read(code_in, fmt=fmt, iostat=ierr) code_out
    if ( ierr .ne. 0 ) then
       ! try and open 'code_in' as a file in dirn
       write(fmt, '(AI0AI0A)') "(A",len_trim(dirn)+1,"A",len_trim(code_in),")"
       write(fname, fmt=fmt) dirn//'/', code_in
       vector_file_lun = vector_file_lun + 1
       open(unit=vector_file_lun, file=fname, iostat=ierr, status='old')
       if ( ierr .ne. 0 ) then
          write(0, *) "Failed to open JPL Ephemeris at ",fname," error: ", ierr
          code_out=0
       else
          code_out=-vector_file_lun
       end if
    end if
    return 

  end subroutine get_code

  subroutine read_sur (dirn, lun_in, point, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine opens and reads in the survey description file.
! Angles are returned in radian.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : October 2004
! Version 3 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     dirn  : Name of directory with survey definition (CH)
!     lun_in: File unit (I4)
!
! OUTPUT
!     point : Pointing structure including all characterisation (pointing)
!     ierr  : Error code
!                0 : nominal run
!               10 : unable to open pointing file
!               20 : error reading record
!               30 : end of file reached
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) dirn
!f2py intent(in) lun_in
!f2py intent(out) point
!f2py intent(out) ierr
    implicit none

    type(t_pointing), intent(out) :: point
    integer, intent(in) :: lun_in
    integer, intent(out) :: ierr
    character(*), intent(in) :: dirn
    type(t_v3d) :: vel
    real (kind=8) :: w, h, ra, dec, r
    integer :: j, nw, lw(nw_max), lun_e, ierr_e, i1, i2, i3, i4
    character(100) :: line, fname
    character(80) :: word(nw_max)
    logical, save :: opened, finished

    data opened /.false./

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
       call create_ears(ra, dec, point%poly)
       do j = nw, 4, -1
          word(j+1) = word(j)
          lw(j+1) = lw(j)
       end do
       nw = nw + 1
    else if (word(1)(1:4) .eq. 'poly') then
       if (nw .lt. 8) goto 2000
       read (word(2), *, err=2000) point%poly%n
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
       do j = 1, point%poly%n
          read (lun_in, *, err=2000, end=3000) point%poly%x(j), point%poly%y(j)
          point%poly%x(j) = point%poly%x(j)*drad
          point%poly%y(j) = point%poly%y(j)*drad
       end do
       call create_poly(ra, dec, point%poly)
!         write (6, *) 'This feature is not implemented yet.'
!         goto 2000
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
       call create_rectangle(w, h, ra, dec, point%poly)
    end if
    call check_polygon(point%poly)
    read (word(5), *, err=2000) point%o_pos(1)%jday
    read (word(6), *, err=2000) point%ff
    ! get the path to the 
    call get_code(word(7), dirn, point%code)
!    read (word(7), *, err=2000) point%code

! USE OF SLALIB: need to get longitude, latitude and elevation of
! observatory. This is given by the sla_OBS routine. One then needs to
! get the LST (see documentation on EXPLANATION AND EXAMPLES:
! Ephemerides).

    point%efnam = word(8)
    call read_file_name (point%efnam, i3, i4, finished, len(point%efnam))

! Open and read in efficiency function
    fname(1:i2-i1+2) = dirn(i1:i2)//'/'
    fname(i2-i1+3:) = point%efnam
    call read_eff (fname, lun_e, point%c, ierr_e)

    if (ierr_e .eq. 10) then
       write (6, *) 'Unable to open '//word(8)
       goto 2000
    else if (ierr_e .eq. 0) then
       goto 1610
    else 
       write (6, *) 'Unknown return code in read_sur.'
       ierr = ierr_e
       return
    end if
1610 continue

! Get rid of unused bins at high magnitude for lookup tables
    do i1 = 1, point%c%nr
       if (point%c%eff_p(i1)%n .gt. 0) then
          j = point%c%eff_p(i1)%n
1700      continue
          if (point%c%eff_p(i1)%e(j) .le. 0.d0) then
             j = j - 1
             goto 1700
          end if
          point%c%eff_p(i1)%n = amin0(j+1, point%c%eff_p(i1)%n)
       end if
    end do

! Computes observatory position at given jday, in ICRF
    call ObsPos (point%code, point%o_pos(1)%jday, point%o_pos(1)%pos, vel, &
         point%o_pos(1)%r, ierr_e)

    if (ierr_e .ne. 0) then
       write (6, *) 'Error while computing observatory''s position. (one)'
       write (6, *) 'ierr = ', ierr_e
       goto 2000
    end if

! The same, 2 hours later
    point%o_pos(2)%jday = point%o_pos(1)%jday + TwoHours
    call ObsPos (point%code, point%o_pos(2)%jday, point%o_pos(2)%pos, vel, &
         point%o_pos(2)%r, ierr_e)
    if (ierr_e .ne. 0) then
       write (6, *) 'Error while computing observatory''s position. (two)'
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

  end subroutine read_sur

  subroutine GetSurvey (survey, lun_s, n_sur, points, sur_mm, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine reads in a survey description.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2004
! Version 2 : October 2004
! Version 3 : January 2006
! Version 4 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     survey: Survey directory name (CH)
!     lun_s : Logical unit for file (I4)
!
! OUTPUT
!     n_sur : Number of pointings read (I4)
!     points: Array of pointings (n*pointing)
!     sur_mm: Limiting magnitude for each survey (n*R8)
!     ierr  : Error code (I4)
!                0 : nominal run
!              100 : Maximum number of objects reached
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) survey
!f2py intent(in) lun_s
!f2py intent(out) n_sur
!f2py intent(out) points
!f2py intent(out) sur_mm
!f2py intent(out) ierr
    implicit none

    type(t_pointing), intent(out) :: points(n_sur_max)
    integer, intent(in) :: lun_s
    integer, intent(out) :: n_sur, ierr
    real (kind=8), intent(out) :: sur_mm(n_sur_max)
    character(*), intent(in) :: survey
    type(t_pointing) :: point
    real (kind=8) :: rate, tmp, mag, eff
    integer :: nr, n, j, i, i1, i2
    logical :: finished

! Open and read in survey definitions
    call read_file_name (survey, i1, i2, finished, len(survey))
    n_sur = 0
200 continue
       call read_sur (survey(i1:i2), lun_s, point, ierr)

       if (ierr .ne. 0) then
          if (ierr .eq. 10) then
             write (6, *) 'Unable to open ',survey(i1:i2),'/pointings.list'
             ierr = -10
          else if (ierr .eq. 20) then
             write (6, *) 'Error reading ',survey(i1:i2),'/pointings.list'
             write (6, *) 'Survey number: ', n_sur
             goto 200
          else if (ierr .eq. 30) then
             goto 300
          else
             write (6, *) 'Unknown return code in read_obj.'
          end if
          return
       end if

       n_sur = n_sur + 1
       points(n_sur) = point
! START comment out in production mode
!       write (18, *) 'Survey number: ', n_sur
!       write (18, *) points(n_sur)%o_pos(1)%jday, points(n_sur)%ff, &
!            points(n_sur)%code
!       write (18, *) points(n_sur)%o_pos(1)%pos%x, &
!            points(n_sur)%o_pos(1)%pos%y, &
!            points(n_sur)%o_pos(1)%pos%z, &
!            points(n_sur)%o_pos(1)%r
!       write (18, *) points(n_sur)%efnam, points(n_sur)%c%nr
! END comment out in production mode
       sur_mm(n_sur) = 0.d0
       nr = point%c%nr
       do j = 1, nr
          n = points(n_sur)%c%eff_p(j)%n
! START comment out in production mode
!          write (18, *) j, points(n_sur)%c%eff_p(j)%min, &
!               points(n_sur)%c%eff_p(j)%max, points(n_sur)%c%eff_p(j)%n, &
!               points(n_sur)%c%eff_p(j)%mag_lim
! END comment out in production mode
          if (n .gt. 0) then
             sur_mm(n_sur) = max(sur_mm(n_sur), points(n_sur)%c%eff_p(j)%b(n))
! START comment out in production mode
!             do i = 1, points(n_sur)%c%eff_p(j)%n
!                write (18, *) j, i, points(n_sur)%c%eff_p(j)%b(i), &
!                     points(n_sur)%c%eff_p(j)%e(i)
!             end do
! END comment out in production mode
          else if ((n .eq. -1) .or. (n .eq. -3)) then
! START comment out in production mode
!             do i = 1, 3
!                write (18, *) j, i, points(n_sur)%c%eff_p(j)%e(i)
!             end do
! END comment out in production mode
          else if ((n .eq. -2) .or. (n .eq. -4)) then
! START comment out in production mode
!             do i = 1, 4
!                write (18, *) j, i, points(n_sur)%c%eff_p(j)%e(i)
!             end do
! END comment out in production mode
          else
             write (6, *) 'Got efficiency function type ', n, j
             write (6, *) 'Should be >0, -1, -2, -3 or -4.'
             ierr = -100
             return
          end if
          if (n .lt. 0) then
             rate = 0.5d0*(points(n_sur)%c%eff_p(j)%min + &
                  points(n_sur)%c%eff_p(j)%max)
! START comment out in production mode
!             write (18, *) 'Rate: ', rate
! END comment out in production mode
             mag = 40.d0
250          continue
                mag = mag - 0.1d0
                eff = eta(point%c%eff_p, nr, mag, rate, tmp)
! START comment out in production mode
!                write (18, *) mag, eff
! END comment out in production mode
                if ((eff .eq. 0.d0) .and. (mag .ge. -0.05d0)) goto 250
260          continue
             sur_mm(n_sur) = max(sur_mm(n_sur), mag+0.1d0)
          end if
! START comment out in production mode
!          write (18, *) n_sur, j, sur_mm(n_sur)
! END comment out in production mode
       end do
! START comment out in production mode
!       write (18, *) n_sur, sur_mm(n_sur)
! END comment out in production mode
       goto 200
300 continue

    ierr = 0
    return
  end subroutine GetSurvey

end module getsur
