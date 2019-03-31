module poly_lib

  use poly_dec

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File Polygon-lib.f
!
! J.-M. Petit  Observatoire de Besançon
! Version 1 :  February 2016
!
! The purpose of this library is to provide polygon-oriented routines.
! The first and most important one is:
!     point_in_polygon (p, poly, n)
! which tels if the point "p" is inside, outside or touching the
! polygon "poly".
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  function point_in_polygon(p, poly)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function checks if point "p" is inside the polygon "poly" using
! the quadrant method.
!
! This moves the given point to the origin and shifts the polygon
! accordingly. Then for each edge of the polygon, calc_walk_summand is
! called. If the sum of all returned values from these calls is +4 or
! -4, the point lies indeed inside the polygon. Otherwise, if a
! PolygonsTouching exception has been caught, the point lies on one of
! he edges of the polygon.
!
! Returns the number of nodes of the polygon, if the point lies inside,
! otherwise 1 if the point lies on the polygon and if not, 0.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     p     : Point (array (x,y)) (2*R8)
!     poly  : Polygon structure (polygon)
!
! OUTPUT
!     point_in_polygon: result of the call (I4)
!                Point inside polygon : n
!                Point on polygon     : 1
!                Point outside polygon: 0
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) p
!f2py intent(in) poly

    implicit none

    integer :: point_in_polygon, i, walk_sum, walk
    integer, parameter :: n_max=100
    type(t_polygon) :: poly
    real (kind=8) :: p(2), moved(2,n_max+1)

! Move point to origin
    do i = 1, poly%n+1
       moved(1,i) = poly%x(i) - p(1)
       moved(2,i) = poly%y(i) - p(2)
    end do

! Computing the running sum
    walk_sum = 0
    do i = 1, poly%n
       walk = calc_walk_summand(moved(1,i), moved(1,i+1))
       if (walk .eq. -100) then
! Point is touching the polygon
          point_in_polygon = 1
          return
       end if
! Point is not on polygon
       walk_sum = walk_sum + walk
    end do

! Final check
    if (abs(walk_sum) .eq. 4) then
       point_in_polygon = poly%n
    else
       point_in_polygon = 0
    end if
    return

  end function point_in_polygon

  function calc_walk_summand(p1, p2)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function calculates the summand along one edge depending on axis
! crossings.
!
! Follows the edge between two points and checks if one or both axes
! are being crossed. If They are crossed in clockwise sense, it returns
! +1 otherwise -1. Going through the origin raises the PolygonsTouching
! exception (returns -100).
!
! Returns one of -2, -1, 0, +1, +2 or raises PolygonsTouching
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     p1    : First point (array (x,y)) of edge (2*R8)
!     p2    : Second point of edge (2*R8)
!
! OUTPUT
!     calc_walk_summand: result of the call (I4)
!                Clockwise crossing        : +1
!                No crossing               : 0
!                Counter-clockwise crossing: -1
!                Diagonal crossing         : +/-2
!                Point on edge             : -100
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) p1
!f2py intent(in) p2

    implicit none

! Indices for better readability
    integer, parameter :: x=1, y=2
    integer :: summand, calc_walk_summand
    real (kind=8) :: p1(2), p2(2), tx, ty, x_y0, y_x0

    summand = 0
! Here, we assume the 2 points are different!
!
! Checking for vertical line
    if (p1(x) .ne. p2(x)) then
       ty = p1(x)/(p1(x) - p2(x))
    else
       ty = p1(y)/(p1(y) - p2(y))
    end if

! Checking for horizontal line
    if (p1(y) .ne. p2(y)) then
       tx = p1(y)/(p1(y) - p2(y))
    else
       tx = ty
    end if

! Compute position of axis intersection
    x_y0 = p1(x) + tx*(p2(x) - p1(x))
    y_x0 = p1(y) + ty*(p2(y) - p1(y))

! Check if crossing x axis
    if ((tx .ge. 0.d0) .and. (tx .lt. 1.d0)) then
! Check if origin on edge
       if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
          calc_walk_summand = -100
          return
       end if
       x_y0 = x_y0*(p2(y) - p1(y))
       if (x_y0 .ne. 0.d0) summand = summand + sign(1.d0, x_y0)
    end if

! Check if crossing y axis
    if ((ty .ge. 0.d0) .and. (ty .lt. 1.d0)) then
! Check if origin on edge
       if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
          calc_walk_summand = -100
          return
       end if
       y_x0 = y_x0*(p1(x) - p2(x))
       if (y_x0 .ne. 0.d0) summand = summand + sign(1.d0, y_x0)
    end if
    calc_walk_summand = summand
    return

  end function calc_walk_summand

  subroutine check_polygon(poly)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! NOTE: user must make sure the input array is at least (2,n+1) long.
! This subroutine will copy first point into index n+1, if last point
! not already same as first. It will also check that there are no two
! consecutive points the same.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     poly  : Polygon structure (polygon)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) poly
    implicit none

    type(t_polygon) :: poly
    integer :: i, j
 
    if ((poly%x(poly%n) .eq. poly%x(1)) .and. (poly%y(poly%n) .eq. poly%y(1))) &
         then
       poly%n = poly%n - 1
    else
       poly%x(poly%n+1) = poly%x(1)
       poly%y(poly%n+1) = poly%y(1)
    end if
    j = 1
1000 continue
    if ((poly%x(j) .eq. poly%x(j+1)) .and. (poly%y(j) .eq. poly%y(j+1))) then
       do i = j+1, poly%n
          poly%x(i) = poly%x(i+1)
          poly%y(i) = poly%y(i+1)
       end do
       poly%n = poly%n - 1
       j = j - 1
    end if
    j = j + 1
    if (j .le. poly%n) goto 1000
    return

  end subroutine check_polygon
end module poly_lib
