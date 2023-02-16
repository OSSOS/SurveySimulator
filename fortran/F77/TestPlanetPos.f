      integer*4 nplanets, ind, i, ierr, j

      real*8 pos_p(3), vel_p(3), pos_b(3), vel_b(3), jday

      jday = 2457023.5d0
      nplanets = 8

      do i = 1, 10
         print *, 'Time =', jday
         print *, 'old:'
         do ind = 1, nplanets
            call PlanetXV(ind, jday, pos_p, vel_p, ierr)
            print *, 'Planet', ind
            print *, (pos_p(j), j=1,3)
         end do
         print *, 'new:'
         do ind = 1, nplanets
            call PlanetXV1(ind, jday, pos_p, vel_p, ierr)
            print *, 'Planet', ind
            print *, (pos_p(j), j=1,3)
         end do
         print *, 'new2:'
         do ind = 1, nplanets
            call PlanetXV2(ind, jday, pos_p, ierr)
            print *, 'Planet', ind
            print *, (pos_p(j), j=1,3)
         end do
         print *, 'Barycenter:'
         print *, 'old:'
         call BaryXV(jday, pos_b, vel_b, ierr)
         print *, (pos_b(j), j=1,3)
         print *, 'new:'
         call BaryXV1(jday, pos_b, vel_b, ierr)
         print *, (pos_b(j), j=1,3)
         print *, ' '
         jday = jday + 100.d0
      end do

      end

      include 'PosVelUtils.f'
