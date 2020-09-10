module poly_dec
  ! define length of array parameters
  integer, parameter :: n_e_max=41

  ! define derived types
  type t_polygon
     integer :: n
     real (kind=8), dimension(n_e_max) :: x, y
  end type t_polygon

end module poly_dec
