program test_read_pointings

  use datadec
  use getsur

  implicit none
  integer, parameter :: lun_s = 13
  character(80) :: survey
  integer :: n_sur, ierr
  real (kind=8) :: sur_mm(n_sur_max)
  type(t_pointing) :: points(n_sur_max)

  survey = './OSSOS'
  call GetSurvey (survey, lun_s, n_sur, points, sur_mm, ierr)

  stop

end program test_read_pointings
