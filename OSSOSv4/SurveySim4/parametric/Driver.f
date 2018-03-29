c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File Driver.f
c
c J.-M. Petit  Observatoire de Besançon
c Version 1 :  May 2013
c
c The purpose of this program is to compare models of the Kuiper Belt to
c reality by actually comparing what the surveys (say CFEPS or OSSOS)
c would have found if the model was a good representation of the real
c world, to the objects that were actually found.
c
c The workflow of the survey simulator driver is as follows:
c
c     Loop (some condition):
c         call GiMeObj(arg_list_1)
c         Check model ended:
c             set exit condition
c         call Detos1(arg_list_2)
c         Check detection and tracking:
c             store results
c
c The GiMeObj routine is in charge of providing a single new object at
c each call. The Detos1 routine determines is the proposed object would
c have been detected by the survey.
c
c The survey simulator expects orbital elements with respect to ecliptic
c reference frame.
c
c Logical unit numbers 7 to 19 are reserved to use by Driver.f and
c SurveySubs.f and should not be used by GiMeObj or any other routine
c that you may add to the driver. Please use logical unit numbers
c starting from 20.
c
c As currently written, the driver includes a file 'GiMeObj.f'
c containing the definition of the model. The correct way to use this
c feature is to have one's GiMeObj routine in a file <whatever.f> and
c create a symobolic link: 
c
c     ln -s <whatever.f> GiMeObj.f
c
c For more information, please refer to ../README.first, ./README.src
c and ./README.surveysubs.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      include 'GiMeObj.f'
      include 'SurveySubs.f'

      implicit none

      integer*4 n_obj_max, screen, keybd, verbose, lun_h, lun_t

      parameter
     $  (n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9)

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)

c color array NEEDS to be length 10 or more!
      real*8 a, e, inc, node, peri, M, h, epoch, m_int, d_ra, d_dec,
     $  r, delta, ra, dec, random, mt, color(10), gb, ph, period, amp,
     $  jday_p, m_rand, eff, rn_iter, eff_lim, h_rand

      integer*4 n_hits, n_track, ierr, seed, flag, isur, ic, n_iter,
     $  n_track_max, nchar

      character distri_file*80, survey_dir*100, 
     $  trk_outfile*80, det_outfile*80,
     $  comments*100, surna*10

      logical keep_going

      lun_h = 10
      lun_t = 11
      keep_going = .true.

c Get arguments
c Seed for random number generator
      read (5, *, err=9999) seed
c Maximum number of detections (>0) or -maximum number of trials (<0)
      read (5, *, err=9999) n_track_max
c Directory containing the characterization files
      read (5, '(a)', err=9999) survey_dir
c File with model parameters
      read (5, '(a)', err=9999) distri_file
c Name for the detected objects outfile
      read (5, '(a)', err=9999) det_outfile
c Name for the tracked objects outfile
      read (5, '(a)', err=9999) trk_outfile

c Open output files and write header
      open (unit=lun_h, file=det_outfile, status='new', err=9500)
      write (lun_h, '(''# Seed: '', i10)') seed
      write (lun_h, '(''#'')')
      write (lun_h, '(a,a,a)')
     $  '#   a      e        i        q        r        M      m_rand',
     $  ' H_rand color flag delta    m_int   H_int eff   RA(H)     DEC',
     $  '    delta_ra delt_dec Surv.  Comments'
      write (lun_h, '(''#'')')
      write (lun_h, '
     $(''# flag: >0: detected; >2: characterized; 0 mod(2): tracked'')')
      write (lun_h, '(''# Survey: name of the block'')')
      write (lun_h, '(''# delta_ra: distance from center of pointing'',
     $  '' [arcsec]'')')
      write (lun_h, '(''# delt_dec: distance from center of pointing'',
     $  '' [arcsec]'')')
      write (lun_h, '(''#'')')

      open (unit=lun_t, file=trk_outfile, status='new', err=9501)
      write (lun_t, '(a,a)')
     $  '#   a      e        i        q        r        M      m_rand',
     $  ' H_rand color Comments'
      write (lun_t, '(''#'')')

c Initialize counters
      n_hits = 0
      n_track = 0
      n_iter = 0
      rn_iter = 0.d0

c Open and read in object distribution
 100  continue
      
c     Main loop: loop on objects
      if (keep_going) then

c        Select object.
         nchar = 0
         call GiMeObj (distri_file, seed, a, e, inc, node, peri,
     $     M, epoch, h, color, gb, ph, period, amp, comments, nchar,
     $     ierr)

         if (ierr .eq. -10) then       !Something wrong with this object, go to next one
            goto 100
         else if (ierr .eq. -20) then     !Something very wrong happend, stop
            write (screen, '(a)') 'GiMeObj returned -20, stopping.'
            goto 2000
         else if (ierr .eq. 100) then     !Reached end of model, prepare to stop
            keep_going = .false.
         end if

c        Count number of iterations; may exceed 2**31, so use a real*8 helper
         n_iter = n_iter + 1
         if (n_iter .gt. 2000000000) then
            rn_iter = rn_iter + dble(n_iter)
            n_iter = 0
            if (rn_iter .ge. 9.d9) goto 2000
         end if

c        Determine if the object would be detected
         call Detos1 (a, e, inc, node, peri, M, epoch, h, color, gb, ph,
     $     period, amp, survey_dir, seed, flag, ra, dec, d_ra, d_dec, r,
     $     delta, m_int, m_rand, eff, isur, mt, jday_p, ic, surna,
     $     h_rand)

c        Check if detected and tracked, and write to output files
         if (flag .gt. 0) then
            n_hits = n_hits + 1
            write (lun_h, 9000) a, e, inc/drad, a*(1.d0-e), r, mt/drad,
     $        m_rand, h_rand, color(ic),
     $        flag, delta, m_int, h, eff, ra/drad/15.,
     $        dec/drad, d_ra/drad*3600./24.,
     $        d_dec/drad*3600./24., surna, comments(1:nchar)
            if ((flag .gt. 2) .and. (mod(flag,2) .eq. 0)) then
               n_track = n_track + 1
               write (lun_t, 9010) a, e, inc/drad, a*(1.d0-e), r,
     $           mt/drad, m_rand, h_rand, color(ic), comments(1:nchar)
            end if
         end if

c Should we continue ?
         if (((n_track_max .gt. 0) .and. (n_track .ge. n_track_max))
     $     .or.
     $     ((n_track_max .lt. 0) .and. (n_iter .ge. -n_track_max))) then
            keep_going = .false.
         end if
         goto 100
c     end of the if( keep_going ) loop    
      end if

 2000 continue
      write (lun_h, '(''#'')')
      write (lun_h, '(''# Total number of objects:   '', f11.0)')
     $  rn_iter + dble(n_iter)
      write (lun_h, '(''# Number of detections:      '', i7)') n_hits
      write (lun_h, '(''# Number of tracked objects: '', i7)') n_track
      close (lun_h)
      close (lun_t)

      stop

 9000 format (f8.3,1x,f6.3,1x,5(f8.3,1x),f6.2,1x,f5.2,1x,i2,2(1x,f8.3),
     $  1x,f6.2,1x,f4.2,1x,f8.5,1x,f8.4,2(1x,f8.5),1x,a6,1x,a)
 9010 format (f8.3,1x,f6.3,1x,5(f8.3,1x),f6.2,1x,f5.2,1x,a)

 9500 continue
      write (screen, *) 'File "', det_outfile, '" already exists. '
      goto 9502

 9501 continue
      write (screen, *) 'File "', trk_outfile, '" already exists. '
      goto 9502

 9502 continue
      write (screen, *) 'Make sure "', det_outfile, '" and "',
     $  trk_outfile, '" do not exist and restart SurveySimulator.'
      stop

 9999 continue
      write (screen, *) 'Usage: SurveySimulator < input'
      write (screen, *)
      write (screen, *)
     $  'This will read the following from the keyboard:'
      write (screen, *) '<seed>'
      write (screen, *) '<n_track_max>'
      write (screen, *) '<survey_dir>'
      write (screen, *) '<distrib_file>'
      write (screen, *) '<det_file>'
      write (screen, *) '<track_file>'
      write (screen, *) 'where:'
      write (screen, *)
     $  '<seed>: integer used as seed for the random number generator'
      write (screen, *)
     $  '<n_track_max>: maximum number of simulated tracked'
     $  //' detection if > 0,'
      write (screen, *)
     $  '              -maximum number of trials if < 0'
      write (screen, *) 
     $  '<survey_dir>: directory of survey charactarization files'
      write (screen, *) 
     $  '<distrib_file>: name of the input file for GiMeObj'
      write (screen, *) 
     $  '<det_file>: output file for detected objects'
      write (screen, *) 
     $  '<track_file>: output file for tracked objects'

      stop

      end
