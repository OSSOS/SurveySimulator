MODULE DE_ephhdr

  USE DE_size

  implicit none
  real (kind=8) :: SS(3), VALS(NMAX), AU, EMRAT
  INTEGER (kind=4) :: IPT(3,15), DENUM, NCON
  data IPT /45*0/
  data DENUM /0/

END MODULE DE_ephhdr
