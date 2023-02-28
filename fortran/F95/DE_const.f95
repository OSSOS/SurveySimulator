MODULE DE_const

  implicit none
  integer (kind=4), PARAMETER :: npl = 10
  real (kind=8) :: masses(15), obli
  integer (kind=4) :: idtarg(15), nctr, ii
  logical :: first

! Sun, Mercury, Venus, Earth, Moon, Mars bary, Jupiter bary, Saturn bary,
! Uranus bary, Neptune Bary, Pluto bary, EMB
  data (idtarg(ii),ii=1,15) /11, 1, 2, 3, 10, 4, 5, 6, 7, 8, 9, 13, 0, 0, 0/
! Center is Solar System barycenter
  data nctr /12/
  data first /.true./

END MODULE DE_const
