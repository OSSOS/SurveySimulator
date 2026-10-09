program test_pointing_geom
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Load CFEPS characterization and check pointing_geom fill factor for
! the first pointing (L3f-smooth.eff, fill=0.80).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  use parameters
  use surveysub

  implicit none

  integer :: n_sur, ierr, nfail
  real (kind=8) :: area, fill
  character(name_len) :: survey
  character(block_len) :: block
  character(key_len) :: key
  character(eff_name_len) :: eff_file
  real (kind=8) :: epoch, mag_lim, rate_mid
  integer :: code

  nfail = 0
  call reset_simulator()
  call survey_load('Surveys/CFEPS', 21, n_sur, ierr)
  if ((ierr .ne. 0) .or. (n_sur .le. 0)) then
     write (6, *) 'FAIL: survey_load ierr=', ierr, ' n_sur=', n_sur
     stop 1
  end if
  write (6, *) 'PASS: loaded n_sur=', n_sur

  call pointing_geom(1, area, fill)
  if (area .le. 0.d0) then
     write (6, *) 'FAIL: area <= 0:', area
     nfail = nfail + 1
  else
     write (6, *) 'PASS: area=', area
  end if
  if (dabs(fill - 0.80d0) .gt. 1.d-6) then
     write (6, *) 'FAIL: fill expected 0.80 got', fill
     nfail = nfail + 1
  else
     write (6, *) 'PASS: fill=', fill
  end if

  call pointing_meta(1, survey, block, key, eff_file, epoch, code, mag_lim, &
       rate_mid)
  if (index(block, 'L3f-smooth') .le. 0) then
     write (6, *) 'FAIL: unexpected block ', block
     nfail = nfail + 1
  else
     write (6, *) 'PASS: block=', block(1:len_trim(block))
  end if
  if (index(key, 'CFEPS/L3f-smooth') .le. 0) then
     write (6, *) 'FAIL: unexpected key ', key
     nfail = nfail + 1
  else
     write (6, *) 'PASS: key=', key(1:len_trim(key))
  end if
  if (index(eff_file, 'L3f-smooth') .le. 0) then
     write (6, *) 'FAIL: unexpected eff_file ', eff_file
     nfail = nfail + 1
  else
     write (6, *) 'PASS: eff_file=', eff_file(1:len_trim(eff_file))
  end if

  if (nfail .eq. 0) then
     write (6, *) 'ALL TESTS PASSED'
     stop 0
  else
     write (6, *) 'FAILURES:', nfail
     stop 1
  end if

end program test_pointing_geom
