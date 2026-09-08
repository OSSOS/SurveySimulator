program test_polygon_area
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Unit tests for polygon_area: a small equatorial rectangle should have
! area approximately width_deg * height_deg.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  use parameters
  use datadec
  use poly_lib
  use getsur

  implicit none

  type(t_polygon) :: poly
  real (kind=8) :: w, h, ra, dec, area, expected, rel
  integer :: nfail

  nfail = 0

  ! 1 deg x 1 deg rectangle at Dec = 0 (half-width/height in radians)
  w = 0.5d0*drad
  h = 0.5d0*drad
  ra = 0.d0
  dec = 0.d0
  call create_rectangle(w, h, ra, dec, poly)
  call check_polygon(poly)
  area = polygon_area(poly)
  expected = 1.d0
  rel = dabs(area - expected)/expected
  if (rel .gt. 1.d-3) then
     write (6, *) 'FAIL: 1x1 deg rectangle area=', area, ' expected~', expected
     nfail = nfail + 1
  else
     write (6, *) 'PASS: 1x1 deg rectangle area=', area
  end if

  ! 2 deg x 3 deg rectangle
  w = 1.0d0*drad
  h = 1.5d0*drad
  call create_rectangle(w, h, ra, dec, poly)
  call check_polygon(poly)
  area = polygon_area(poly)
  expected = 6.d0
  rel = dabs(area - expected)/expected
  if (rel .gt. 1.d-3) then
     write (6, *) 'FAIL: 2x3 deg rectangle area=', area, ' expected~', expected
     nfail = nfail + 1
  else
     write (6, *) 'PASS: 2x3 deg rectangle area=', area
  end if

  ! Degenerate polygon
  poly%n = 2
  area = polygon_area(poly)
  if (area .ne. 0.d0) then
     write (6, *) 'FAIL: n<3 should give 0, got', area
     nfail = nfail + 1
  else
     write (6, *) 'PASS: n<3 returns 0'
  end if

  if (nfail .eq. 0) then
     write (6, *) 'ALL TESTS PASSED'
     stop 0
  else
     write (6, *) 'FAILURES: ', nfail
     stop 1
  end if

end program test_polygon_area
