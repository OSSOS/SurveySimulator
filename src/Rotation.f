      subroutine equ_ecl(epsilon, poseq, posecl)
c******************************************************
c
c    Transformation of a vector from the equatorial frame to the ecliptic frame
c
c    ANGLE IN RADIAN !!!
c
c  Author  : F. Mignard OCA/CERGA
c
c
c*** INPUT
c       epsilon : obliquity
c       poseq   : input vector in equatorial frame
c
c*** OUTPUT
c       posecl  : output vector in ecliptic frame
c
c******************************************************
      implicit none

      real*8
     $  Pi, degrad, phi, poseq(3), posecl(3), epsilon, phir

c Chapront et al. 2002 gamma to O_icrs in arcsec
      parameter
     $  (Pi = 3.141592653589793238d0, phi = 0.05542d0,
     $  degrad = Pi/180.d0)

      phir    = phi/3600d0*degrad
 
      call RotZ(-phir, poseq, posecl)

      call RotX(epsilon, posecl, posecl)

      return
      end

      subroutine ecl_equ(epsilon, posecl, poseq)
c******************************************************
c
c    Transformation of a vector from the ecliptic frame to the equatorial frame
c
c    ANGLE IN RADIAN !!!
c
c  Author  : F. Mignard OCA/CERGA
c
c
c*** INPUT
c       epsilon : obliquity
c       poseq   : input vector in ecliptic frame
c
c*** OUTPUT
c       posecl  : output vector in equatorial frame
c
c******************************************************
      implicit none

      real*8
     $  Pi, degrad, phi, poseq(3), posecl(3), epsilon, phir

c Chapront et al. 2002 gamma to O_icrs in arcsec
      parameter
     $  (Pi = 3.141592653589793238d0, phi = 0.05542d0,
     $  degrad = Pi/180.d0)

      phir    = phi/3600d0*degrad

      call RotX(-epsilon, posecl, poseq)

      call RotZ(phir, poseq, poseq)

      return
      end

      subroutine RotX (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the X axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(1) = pos(1)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(2) = pos(2)*ca + pos(3)*sa
      posout(3) = pos(3)*ca - pos(2)*sa

      return
      end

      subroutine RotY (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the Y axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(2) = pos(2)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(3) = pos(3)*ca + pos(1)*sa
      posout(1) = pos(1)*ca - pos(3)*sa

      return
      end

      subroutine RotZ (alpha, posin, posout)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes position in a frame rotated by angle alpha (rd)
c about the Z axis.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : June 2005
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     alpha : Rotation angle, in radian (R8)
c     posin : Position in original frame (R8)
c
c OUTPUT
c     posout: Position in new frame (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer
     $  i

      real*8
     $  alpha, posin(*), posout(*), pos(3), ca, sa

      do i = 1, 3
         pos(i) = posin(i)
      end do
      posout(3) = pos(3)
      ca = dcos(alpha)
      sa = dsin(alpha)
      posout(1) = pos(1)*ca + pos(2)*sa
      posout(2) = pos(2)*ca - pos(1)*sa

      return
      end

      subroutine equat_ecl(ieqec,v_in,v_out,ierr)
!
!     Transformation of a vector v_in(3) from :
!        equator to ecliptic  : irot = +1
!        ecliptic to equator  : irot = -1
!        at J2000.
!     The equator is assumed to be the ICRF frame
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call equat_ecl(ieqec,vv,vv,ierr)
!
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from equator to ecliptic ; -1 from ecliptic to equator.
!     v_in     :  input vector v_in(3) in the initial frame
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
Cf2py intent(in) ieqec
Cf2py intent(in) v_in
Cf2py intent(out) v_out
Cf2py intent(out) ierr

      implicit none

      real*8
     $  epsilon, v_in(3), v_out(3), ww(3), coseps, sineps, Pi, secrad

c obliquity at J20000 in arcsec

      parameter
     $  (epsilon = 84381.41d0, Pi = 3.141592653589793238d0,
     $  secrad = Pi/180.d0/3600.d0)

      integer*4
     $  ieqec, ierr

      coseps = dcos(epsilon*secrad)
      sineps = dsin(epsilon*secrad)

c regular exit
      ierr   = 0

c to allow a call like ::  call equat_ecl(ieqec,vv,vv,ierr)
      ww(1)  = v_in(1)
      ww(2)  = v_in(2)
      ww(3)  = v_in(3)

      if (ieqec .eq. 1) then
        v_out(1) =   ww(1)
        v_out(2) =   coseps*ww(2) + sineps*ww(3)
        v_out(3) = - sineps*ww(2) + coseps*ww(3)
      else if (ieqec .eq. -1) then
        v_out(1) =   ww(1)
        v_out(2) =   coseps*ww(2) - sineps*ww(3)
        v_out(3) =   sineps*ww(2) + coseps*ww(3)
      else
c anomalous exit ieqec not allowed
        ierr     =   100
      end if

      return
      end

      subroutine invar_ecl(ieqec,v_in,v_out,ierr)
!
!     Transformation of a vector v_in(3) from :
!        invariable plane to ecliptic  : irot = +1
!        ecliptic to invariable plane  : irot = -1
!        at J2000.
!     The invariable plane is given with respect to the ecliptic plane
!     (J2000) by Burkhardt, AA, 1982:
!     inclination of invariable plane:   1° 35' 13.86" =   5713.86"
!     direction of ascending node:     107° 36' 30.8"  = 387390.8"
!
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call equat_ecl(ieqec,vv,vv,ierr)
!
!     Adapted from:
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!     J-M. Petit  UBC/CNRS
!     Version 1 : September 2007
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from equator to ecliptic ; -1 from ecliptic to equator.
!     v_in     :  input vector v_in(3) in the initial frame
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
Cf2py intent(in) ieqec
Cf2py intent(in) v_in
Cf2py intent(out) v_out
Cf2py intent(out) ierr

      implicit none

      real*8
     $  epsilon, v_in(3), v_out(3), ww(3), coseps, sineps, Pi, secrad,
     $  omega, cosom, sinom

c obliquity at J20000 in arcsec

      parameter
     $  (epsilon = 5713.86d0, omega = 387390.8d0,
     $  Pi = 3.141592653589793238d0, secrad = Pi/180.d0/3600.d0)

      integer*4
     $  ieqec, ierr

      coseps = dcos(epsilon*secrad)
      sineps = dsin(epsilon*secrad)
      cosom = dcos(omega*secrad)
      sinom = dsin(omega*secrad)

c regular exit
      ierr   = 0

c to allow a call like ::  call equat_ecl(ieqec,vv,vv,ierr)
      ww(1)  = v_in(1)
      ww(2)  = v_in(2)
      ww(3)  = v_in(3)

      if (ieqec .eq. 1) then
        v_out(1) =   cosom*ww(1) - sinom*(coseps*ww(2) - sineps*ww(3))
        v_out(2) =   sinom*ww(1) + cosom*(coseps*ww(2) - sineps*ww(3))
        v_out(3) =   sineps*ww(2) + coseps*ww(3)
      else if (ieqec .eq. -1) then
        v_out(1) =   cosom*ww(1) + sinom*ww(2)
        v_out(2) =   coseps*(-sinom*ww(1) + cosom*ww(2)) + sineps*ww(3)
        v_out(3) = - sineps*(-sinom*ww(1) + cosom*ww(2)) + coseps*ww(3)
      else
c anomalous exit ieqec not allowed
        ierr     =   100
      end if

      return
      end

      subroutine invar_ecl_osc(ieqec, ai, ei, ii, noi, pei, mi,
     $  ao, eo, io, noo, peo, mo, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine convert osculating elements back and forth between
c invariable plane and ecliptic plane.
c Uses invar_ecl to do the work.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
 
      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  posi(3), poso(3), veli(3), velo(3), mu, Pi, TwoPi, drad,
     $  aid, eid, iid, noid, peid, mid

      integer*4
     $  ieqec, ierr

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.d0*Pi, drad = Pi/180.d0,
     $  mu = TwoPi**2)

      aid = ai
      eid = ei
      iid = ii
      noid = noi
      peid = pei
      mid = mi
      call coord_cart (mu, aid, eid, iid, noid, peid, mid, posi(1),
     $  posi(2), posi(3), veli(1), veli(2), veli(3))
      call invar_ecl (ieqec, posi, poso, ierr)
      call invar_ecl (ieqec, veli, velo, ierr)
      call osc_el (mu, poso(1), poso(2), poso(3), velo(1), velo(2),
     $  velo(3), ao, eo, io, noo, peo, mo)
      call ztopi (io)
      if (io .gt. Pi) io = io - Pi
      call ztopi (noo)
      call ztopi (peo)
      call ztopi (mo)

      return
      end

      subroutine invar_ecl_inc_node(ieqec, ii, noi, io, noo, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine convert inclination and longitude of node elements back
c and forth between invariable plane and ecliptic plane.
c Uses invar_ecl_osc to do the work.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
 
      implicit none

      real*8
     $  ai, ei, ii, noi, pei, mi, ao, eo, io, noo, peo, mo,
     $  iid, noid

      integer*4
     $  ieqec, ierr

      data
     $  ai /10.d0/,
     $  ei /0.2d0/,
     $  pei /0.d0/,
     $  mi /0.d0/

      iid = ii
      noid = noi
      call invar_ecl_osc (ieqec, ai, ei, iid, noid, pei, mi, ao, eo, io,
     $  noo, peo, mo, ierr)

      return
      end

      subroutine ztopi (var)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This subroutine resets var to be between 0 and 2Pi.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      real*8
     $  Pi, TwoPi, var

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

 1000 continue
      if (var .gt. TwoPi) then
         var = var - TwoPi
         goto 1000
      end if
 1100 continue
      if (var .lt. 0.d0) then
         var = var + TwoPi
         goto 1100
      end if

      return
      end
