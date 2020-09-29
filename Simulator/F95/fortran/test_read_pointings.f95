program test_read_pointings

  use datadec
  use getsur

  implicit none
  integer, parameter :: lun_s = 13
  character(80) :: survey
  integer :: n_sur, ierr, i
  real (kind=8) :: sur_mm(n_sur_max)
  type(t_pointing) :: points(n_sur_max)

  survey = './SS_Input_Formats'
  call GetSurvey (survey, lun_s, n_sur, points, sur_mm, ierr)
  print *, 'Number of pointings read ', n_sur
  do i = 1, n_sur
     print *, i
     call pr_point(points(i))
  end do

  stop

contains

  subroutine pr_point(point)

    implicit none
    type(t_pointing) :: point

    print *, 'Filling factor: ', point%ff
    print *, 'Observatory code: ', point%code
    print *, 'Efficiency filename: ', point%efnam
    print *, 'Julian Day (1): ', point%o_pos(1)%jday
    print *, 'Sun distance (1): ', point%o_pos(1)%r
    print *, 'Obs position (1): ', point%o_pos(1)%pos
    print *, 'Julian Day (2): ', point%o_pos(2)%jday
    print *, 'Sun distance (2): ', point%o_pos(2)%r
    print *, 'Obs position (2): ', point%o_pos(2)%pos
    print *, 'Number of points in polygon: ', point%poly%n
    print *, 'Filter index: ', point%c%f
    print *, 'Number of rate ranges: ', point%c%nr
  end subroutine pr_point

end program test_read_pointings
