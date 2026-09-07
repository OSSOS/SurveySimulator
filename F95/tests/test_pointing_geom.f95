program test_pointing_geom
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Load CFEPS characterization and check pointing_geom fill factor for
! the first pointing (L3f-smooth.eff, fill=0.80).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  use surveysub

  implicit none

  integer :: n_sur, ierr, nfail
  real (kind=8) :: area, fill
  character(80) :: efnam
  real (kind=8) :: epoch, mag_lim, rate_mid
  integer :: code

  nfail = 0
  call reset_simulator()
  call survey_load('../../src/ossssim/Characterizations/CFEPS', 21, n_sur, ierr)
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

  call pointing_meta(1, efnam, epoch, code, mag_lim, rate_mid)
  if (index(efnam, 'L3f-smooth') .le. 0) then
     write (6, *) 'FAIL: unexpected efnam ', efnam
     nfail = nfail + 1
  else
     write (6, *) 'PASS: efnam=', efnam
  end if

  if (nfail .eq. 0) then
     write (6, *) 'ALL TESTS PASSED'
     stop 0
  else
     write (6, *) 'FAILURES:', nfail
     stop 1
  end if

end program test_pointing_geom
