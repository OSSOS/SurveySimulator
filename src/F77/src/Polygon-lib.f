c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File Polygon-lib.f
c
c J.-M. Petit  Observatoire de Besançon
c Version 1 :  February 2016
c
c The purpose of this library is to provide polygon-oriented routines.
c The first and most important one is:
c     point_in_polygon (p, poly, n)
c which tels if the point "p" is inside, outside or touching the
c polygon "poly".
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      integer*4 function point_in_polygon(p, poly, n)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function checks if point "p" is inside the polygon "poly" using
c the quadrant method.
c
c This moves the given point to the origin and shifts the polygon
c accordingly. Then for each edge of the polygon, calc_walk_summand is
c called. If the sum of all returned values from these calls is +4 or
c -4, the point lies indeed inside the polygon. Otherwise, if a
c PolygonsTouching exception has been caught, the point lies on one of
c he edges of the polygon.
c
c Returns the number of nodes of the polygon, if the point lies inside,
c otherwise 1 if the point lies on the polygon and if not, 0.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     p     : Point (array (x,y)) (2*R8)
c     poly  : Array of points ((2,n+1)*R8)
c     n     : Number of vertices in the polygon (I4)
c
c OUTPUT
c     point_in_polygon: result of the call (I4)
c                Point inside polygon : n
c                Point on polygon     : 1
c                Point outside polygon: 0
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) p
Cf2py intent(in) poly
Cf2py intent(in) n

      implicit none

      integer*4
     $  n, n_max, i, walk_sum, walk, calc_walk_summand

      parameter
     $  (n_max = 100)

      real*8
     $  p(2), poly(2,n+1), moved(2,n_max+1)

      external
     $  calc_walk_summand

c Move point to origin
      do i = 1, n+1
         moved(1,i) = poly(1,i) - p(1)
         moved(2,i) = poly(2,i) - p(2)
      end do

c Computing the running sum
      walk_sum = 0
      do i = 1, n
         walk = calc_walk_summand(moved(1,i), moved(1,i+1))
         if (walk .eq. -100) then
c Point is touching the polygon
            point_in_polygon = 1
            return
         end if
c Point is not on polygon
         walk_sum = walk_sum + walk
      end do

c Final check
      if (abs(walk_sum) .eq. 4) then
         point_in_polygon = n
      else
         point_in_polygon = 0
      end if
      return

      end

      integer*4 function calc_walk_summand(p1, p2)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function calculates the summand along one edge depending on axis
c crossings.
c
c Follows the edge between two points and checks if one or both axes
c are being crossed. If They are crossed in clockwise sense, it returns
c +1 otherwise -1. Going through the origin raises the PolygonsTouching
c exception (returns -100).
c
c Returns one of -2, -1, 0, +1, +2 or raises PolygonsTouching
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     p1    : First point (array (x,y)) of edge (2*R8)
c     p2    : Second point of edge (2*R8)
c
c OUTPUT
c     calc_walk_summand: result of the call (I4)
c                Clockwise crossing        : +1
c                No crossing               : 0
c                Counter-clockwise crossing: -1
c                Diagonal crossing         : +/-2
c                Point on edge             : -100
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) p1
Cf2py intent(in) p2

      implicit none

c Indices for better readability
      integer*4
     $  x, y, summand

      parameter
     $  (x = 1, y = 2)

      real*8
     $  p1(2), p2(2), tx, ty, x_y0, y_x0

      summand = 0
c Here, we assume the 2 points are different!
c
c Checking for vertical line
      if (p1(x) .ne. p2(x)) then
         ty = p1(x)/(p1(x) - p2(x))
      else
         ty = p1(y)/(p1(y) - p2(y))
      end if

c Checking for horizontal line
      if (p1(y) .ne. p2(y)) then
         tx = p1(y)/(p1(y) - p2(y))
      else
         tx = ty
      end if

c Compute position of axis intersection
      x_y0 = p1(x) + tx*(p2(x) - p1(x))
      y_x0 = p1(y) + ty*(p2(y) - p1(y))

c Check if crossing x axis
      if ((tx .ge. 0.d0) .and. (tx .lt. 1.d0)) then
c Check if origin on edge
         if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
            calc_walk_summand = -100
            return
         end if
         x_y0 = x_y0*(p2(y) - p1(y))
         if (x_y0 .ne. 0.d0) summand = summand + sign(1.d0, x_y0)
      end if

c Check if crossing y axis
      if ((ty .ge. 0.d0) .and. (ty .lt. 1.d0)) then
c Check if origin on edge
         if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
            calc_walk_summand = -100
            return
         end if
         y_x0 = y_x0*(p1(x) - p2(x))
         if (y_x0 .ne. 0.d0) summand = summand + sign(1.d0, y_x0)
      end if
      calc_walk_summand = summand
      return

      end

      subroutine check_polygon(poly, n)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c NOTE: user must make sure the input array is at least (2,n+1) long.
c This subroutine will copy first point into index n+1, if last point
c not already same as first. It will also check that there are no two
c consecutive points the same.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : February 2016
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     poly  : Array of points ((2,n+1)*R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      integer*4
     $  n, i, j

      real*8
     $  poly(2,*)

      if ((poly(1,n) .eq. poly(1,1)) .and. (poly(2,n) .eq. poly(2,1)))
     $  then
         n = n - 1
      else
         poly(1,n+1) = poly(1,1)
         poly(2,n+1) = poly(2,1)
      end if
      j = 1
 1000 continue
      if ((poly(1,j) .eq. poly(1,j+1))
     $  .and. (poly(2,j) .eq. poly(2,j+1))) then
         do i = j+1, n
            poly(1,i) = poly(1,i+1)
            poly(2,i) = poly(2,i+1)
         end do
         n = n - 1
         j = j - 1
      end if
      j = j + 1
      if (j .le. n) goto 1000
      return

      end
