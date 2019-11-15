    program check_JMP
!
!   comparison of various sources of solar system ephemeris
!
   use jmp_planetlib
!
   implicit none

   integer, parameter  ::  nstep  = 100! number of steps
   integer, parameter  ::  iout   = 10
   integer, parameter  ::  iframe = 0   ! 0 : equatorial, 1 : ecliptic
   
   

   real(kind=dp), parameter   ::  step  = 400d0   ! step in days
   real(kind=dp), parameter   ::  tbeg  =  xjd2000 ! starting date

 

   integer                      ::  i, iplan

   real(kind=dp)                :: datjd
   real(kind=dp), dimension(3)  :: rr, vv
   real(kind=dp)                :: xlong, xlat, dist 
 
   open(unit = iout, file = 'Results/res_test_JMP.txt', status = 'unknown')

!   do iplan = 0,8  ! explore sun and planets  (1,8 for heliocentric)
   do iplan = 1,8  ! explore planets  (0,8 - sun and planets - for barycentric)
      write(iout, *)
      write(iout, *)  ssname(iplan)
      write(iout, *) '  day-J2010           xx                 yy                   zz                 vx                   vy                  vz   '        
      do i = -nstep, nstep  ! time sampling
!
         datjd = tbeg + i*step 
! Barycentric
!         call ephem_planet_simon(datjd, iplan, iframe, rr, vv)
!         write(iout,'(f12.4,  6e20.10)')  (datjd-tbeg),  rr, vv
! Heliocentric 
         call planet_helio_simon(0, iplan, datjd, xlong, xlat, dist, rr, vv)
         write(iout,'(f12.4,  6e20.10)')  (datjd-tbeg),  rr*xau, vv*audtokms 
 
      enddo
   enddo
   
   print * , ' Program terminated'


    end program check_JMP


   
