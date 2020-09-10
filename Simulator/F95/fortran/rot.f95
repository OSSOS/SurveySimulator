module rot

  use datadec
  use elemutils

contains
  subroutine equ_ecl(epsilon, poseq, posecl)
!******************************************************
!
!    Transformation of a vector from the equatorial frame to the ecliptic frame
!
!    ANGLE IN RADIAN !!!
!
!  Author  : F. Mignard OCA/CERGA
!
!
!*** INPUT
!       epsilon : obliquity
!       poseq   : input vector in equatorial frame
!
!*** OUTPUT
!       posecl  : output vector in ecliptic frame
!
!******************************************************
!f2py intent(in) epsilon
!f2py intent(in) poseq
!f2py intent(out) posecl
    implicit none

    real (kind=8), intent(in) :: poseq(3), epsilon
    real (kind=8), intent(out) :: posecl(3)
! Chapront et al. 2002 gamma to O_icrs in arcsec
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, phi = 0.05542d0, &
         degrad = Pi/180.d0
    real (kind=8) :: phir

    phir    = phi/3600d0*degrad
 
    call RotZ(-phir, poseq, posecl)

    call RotX(epsilon, posecl, posecl)

    return
  end subroutine equ_ecl

  subroutine ecl_equ(epsilon, posecl, poseq)
!******************************************************
!
!    Transformation of a vector from the ecliptic frame to the equatorial frame
!
!    ANGLE IN RADIAN !!!
!
!  Author  : F. Mignard OCA/CERGA
!
!
!*** INPUT
!       epsilon : obliquity
!       posecl  : input vector in ecliptic frame
!
!*** OUTPUT
!       poseq   : output vector in equatorial frame
!
!******************************************************
!f2py intent(in) epsilon
!f2py intent(in) posecl
!f2py intent(out) poseq
      implicit none

    real (kind=8), intent(in) :: posecl(3), epsilon
    real (kind=8), intent(out) :: poseq(3)
! Chapront et al. 2002 gamma to O_icrs in arcsec
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, phi = 0.05542d0, &
         degrad = Pi/180.d0
    real (kind=8) :: phir

    phir    = phi/3600d0*degrad

    call RotX(-epsilon, posecl, poseq)

    call RotZ(phir, poseq, poseq)

    return
  end subroutine ecl_equ

  subroutine RotX (alpha, posin, posout)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes position in a frame rotated by angle alpha (rd)
! about the X axis.
! This is the same as giving the coordinates of a vector that have been
! rotating by angle -alpha around the axis.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2005
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     alpha : Rotation angle, in radian (R8)
!     posin : Position in original frame (R8)
!
! OUTPUT
!     posout: Position in new frame (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) alpha
!f2py intent(in) posin
!f2py intent(out) posout
    implicit none

    real (kind=8), intent(in) :: alpha, posin(3)
    real (kind=8), intent(out) :: posout(3)
    real (kind=8) :: pos(3), ca, sa

    pos = posin

    posout(1) = pos(1)
    ca = dcos(alpha)
    sa = dsin(alpha)
    posout(2) = pos(2)*ca + pos(3)*sa
    posout(3) = pos(3)*ca - pos(2)*sa

    return
  end subroutine RotX

  subroutine RotY (alpha, posin, posout)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes position in a frame rotated by angle alpha (rd)
! about the Y axis.
! This is the same as giving the coordinates of a vector that have been
! rotating by angle -alpha around the axis.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2005
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     alpha : Rotation angle, in radian (R8)
!     posin : Position in original frame (R8)
!
! OUTPUT
!     posout: Position in new frame (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) alpha
!f2py intent(in) posin
!f2py intent(out) posout
    implicit none

    real (kind=8), intent(in) :: alpha, posin(3)
    real (kind=8), intent(out) :: posout(3)
    real (kind=8) :: pos(3), ca, sa

    pos = posin
    posout(2) = pos(2)
    ca = dcos(alpha)
    sa = dsin(alpha)
    posout(3) = pos(3)*ca + pos(1)*sa
    posout(1) = pos(1)*ca - pos(3)*sa

    return
  end subroutine RotY

  subroutine RotZ (alpha, posin, posout)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes position in a frame rotated by angle alpha (rd)
! about the Z axis.
! This is the same as giving the coordinates of a vector that have been
! rotating by angle -alpha around the axis.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2005
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     alpha : Rotation angle, in radian (R8)
!     posin : Position in original frame (R8)
!
! OUTPUT
!     posout: Position in new frame (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) alpha
!f2py intent(in) posin
!f2py intent(out) posout
    implicit none

    real (kind=8), intent(in) :: alpha, posin(3)
    real (kind=8), intent(out) :: posout(3)
    real (kind=8) :: pos(3), ca, sa

    pos = posin
    posout(3) = pos(3)
    ca = dcos(alpha)
    sa = dsin(alpha)
    posout(1) = pos(1)*ca + pos(2)*sa
    posout(2) = pos(2)*ca - pos(1)*sa

    return
  end subroutine RotZ

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
!f2py intent(in) ieqec
!f2py intent(in) v_in
!f2py intent(out) v_out
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(in) :: v_in
    type(t_v3d), intent(out) :: v_out
    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
! obliquity at J20000 in arcsec
    real (kind=8), parameter :: epsilon = 84381.41d0, &
         Pi = 3.141592653589793238d0, secrad = Pi/180.d0/3600.d0
    real (kind=8) :: ww(3), coseps, sineps

    coseps = dcos(epsilon*secrad)
    sineps = dsin(epsilon*secrad)

! regular exit
    ierr   = 0

! to allow a call like ::  call equat_ecl(ieqec,vv,vv,ierr)
    ww(1) = v_in%x
    ww(2) = v_in%y
    ww(3) = v_in%z

    if (ieqec .eq. 1) then
       v_out%x =   ww(1)
       v_out%y =   coseps*ww(2) + sineps*ww(3)
       v_out%z = - sineps*ww(2) + coseps*ww(3)
    else if (ieqec .eq. -1) then
       v_out%x =   ww(1)
       v_out%y =   coseps*ww(2) - sineps*ww(3)
       v_out%z =   sineps*ww(2) + coseps*ww(3)
    else
! anomalous exit ieqec not allowed
       ierr     =   100
    end if

    return
  end subroutine equat_ecl

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
!     One can call invar_ecl(ieqec,vv,vv,ierr)
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
!     Ieqec    : +1 ==> from invariable to ecliptic ; -1 from ecliptic to invariable.
!     v_in     :  input vector v_in(3) in the initial frame
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
!f2py intent(in) ieqec
!f2py intent(in) v_in
!f2py intent(out) v_out
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(in) :: v_in
    type(t_v3d), intent(out) :: v_out
    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
! obliquity at J20000 in arcsec
    real (kind=8), parameter :: epsilon = 5713.86d0, omega = 387390.8d0, &
         Pi = 3.141592653589793238d0, secrad = Pi/180.d0/3600.d0
    real (kind=8) :: ww(3), coseps, sineps, cosom, sinom

    coseps = dcos(epsilon*secrad)
    sineps = dsin(epsilon*secrad)
    cosom = dcos(omega*secrad)
    sinom = dsin(omega*secrad)

! regular exit
    ierr   = 0

! to allow a call like ::  call invar_ecl(ieqec,vv,vv,ierr)
    ww(1) = v_in%x
    ww(2) = v_in%y
    ww(3) = v_in%z

    if (ieqec .eq. 1) then
       v_out%x =   cosom*ww(1) - sinom*(coseps*ww(2) - sineps*ww(3))
       v_out%y =   sinom*ww(1) + cosom*(coseps*ww(2) - sineps*ww(3))
       v_out%z =   sineps*ww(2) + coseps*ww(3)
    else if (ieqec .eq. -1) then
       v_out%x =   cosom*ww(1) + sinom*ww(2)
       v_out%y =   coseps*(-sinom*ww(1) + cosom*ww(2)) + sineps*ww(3)
       v_out%z = - sineps*(-sinom*ww(1) + cosom*ww(2)) + coseps*ww(3)
    else
! anomalous exit ieqec not allowed
       ierr     =   100
    end if

    return
  end subroutine invar_ecl

  subroutine invar_ecl_osc(ieqec, o_mi, o_mo, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine convert osculating elements back and forth between
! invariable plane and ecliptic plane.
! Uses invar_ecl to do the work.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ieqec
!f2py intent(in) o_mi
!f2py intent(out) o_mo
!f2py intent(out) ierr
 
    implicit none

    type(t_orb_m), intent(in) :: o_mi
    type(t_orb_m), intent(out) :: o_mo
    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi, &
         drad = Pi/180.d0, mu = TwoPi**2
    type(t_orb_m) :: o_md
    type(t_v3d) :: posi, poso, veli, velo
    real (kind=8) ::  aid, eid, iid, noid, peid, mid

    o_md%a = o_mi%a
    o_md%e = o_mi%e
    o_md%inc = o_mi%inc
    o_md%node = o_mi%node
    o_md%peri = o_mi%peri
    o_md%m = o_mi%m
    call coord_cart (mu, o_md, posi, veli)
    call invar_ecl (ieqec, posi, poso, ierr)
    call invar_ecl (ieqec, veli, velo, ierr)
    call osc_el (mu, poso, velo, o_mo)
    call ztopi (o_mo%inc)
    if (o_mo%inc .gt. Pi) o_mo%inc = o_mo%inc - Pi
    call ztopi (o_mo%node)
    call ztopi (o_mo%peri)
    call ztopi (o_mo%m)

    return
  end subroutine invar_ecl_osc

  subroutine invar_ecl_inc_node(ieqec, ii, noi, io, noo, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine convert inclination and longitude of node elements back
! and forth between invariable plane and ecliptic plane.
! Uses invar_ecl_osc to do the work.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ieqec
!f2py intent(in) ii
!f2py intent(in) noi
!f2py intent(out) io
!f2py intent(out) noo
!f2py intent(out) ierr
     implicit none

    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: ii, noi
    real (kind=8), intent(out) :: io, noo
    type(t_orb_m) :: o_mi, o_mo

    o_mi%a = 10.d0
    o_mi%e  =0.2d0
    o_mi%peri = 0.d0
    o_mi%m = 0.d0
    o_mi%inc = ii
    o_mi%node = noi
    call invar_ecl_osc (ieqec, o_mi, o_mo, ierr)
    io = o_mo%inc
    noo = o_mo%node

    return
  end subroutine invar_ecl_inc_node

  subroutine ref_ecl(ieqec, v_in, v_out, eps, om, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine converts vectors back and forth between
! given reference frame and ecliptic plane.
! Reference frame is given by Omega (om), the node longitude (rotation
! around Z axis of ecliptic) and Epsilon (eps), the inclination of the
! reference frame (rotation around node axis of ecliptic).
!
!    ANGLE IN RADIAN !!!
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!     Transformation of a vector v_in(3) from :
!        given reference frame to ecliptic  : irot = +1
!        ecliptic to given reference frame  : irot = -1
!        at J2000.
!
!     The ecliptic is the so called 'conventional ecliptic'
!     going through the origin of the ICRF  with the obliquity
!     epsilon = 23°26'21".410 = 84381".41
!     It differs by ~ 50 mas from the inertial mean ecliptic(s) of J2000.
!
!     One can call ref_ecl(ieqec,vv,vv,ierr)
!
!     Adapted from:
!     F. Mignard  OCA/CERGA
!     Version 1 : April 2003
!
!     J-M. Petit  UBC/CNRS
!     Version 1 : June 2015
!
!***************************************************************
! INPUT
!     Ieqec    : +1 ==> from reference frame to ecliptic ; 
!                -1 ==> from ecliptic to reference frame
!     v_in     :  input vector v_in(3) in the initial frame
!     eps      :  Epsilon, inclination of ref frame [rad]
!     om       :  Omega, longitude of node [rad]
!
! OUTPUT
!     v_out    :  output vector in the final frame v_out(3)
!     ierr     :  ierr = 0   :: normal exit
!                 ierr = 100 :: error ieqec neither 1 or -1 . 
!
!***************************************************************
!f2py intent(in) ieqec
!f2py intent(in) v_in
!f2py intent(in) eps
!f2py intent(in) om
!f2py intent(out) v_out
!f2py intent(out) ierr
    implicit none

    type(t_v3d), intent(in) :: v_in
    type(t_v3d), intent(out) :: v_out
    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: eps, om
! obliquity at J20000 in arcsec
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         secrad = Pi/180.d0/3600.d0
    real (kind=8) :: ww(3), coseps, sineps, cosom, sinom

    coseps = dcos(eps)
    sineps = dsin(eps)
    cosom = dcos(om)
    sinom = dsin(om)

! regular exit
    ierr   = 0

! to allow a call like ::  call ref_ecl(ieqec,vv,vv,ierr)
    ww(1) = v_in%x
    ww(2) = v_in%y
    ww(3) = v_in%z

    if (ieqec .eq. 1) then
       v_out%x =   cosom*ww(1) - sinom*(coseps*ww(2) - sineps*ww(3))
       v_out%y =   sinom*ww(1) + cosom*(coseps*ww(2) - sineps*ww(3))
       v_out%z =   sineps*ww(2) + coseps*ww(3)
    else if (ieqec .eq. -1) then
       v_out%x =   cosom*ww(1) + sinom*ww(2)
       v_out%y =   coseps*(-sinom*ww(1) + cosom*ww(2)) + sineps*ww(3)
       v_out%z = - sineps*(-sinom*ww(1) + cosom*ww(2)) + coseps*ww(3)
    else
! anomalous exit ieqec not allowed
       ierr     =   100
    end if

    return
  end subroutine ref_ecl

  subroutine ref_ecl_osc(ieqec, o_mi, o_mo, eps, om, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine converts osculating elements back and forth between
! given reference frame and ecliptic plane.
! Uses ref_ecl to do the work.
! Reference frame is given by Omega (om), the node longitude (rotation
! around Z axis of ecliptic) and Epsilon (eps), the inclination of the
! reference frame (rotation around node axis of ecliptic).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ieqec
!f2py intent(in) o_mi
!f2py intent(in) eps
!f2py intent(in) om
!f2py intent(out) o_mo
!f2py intent(out) ierr
    implicit none

    type(t_orb_m), intent(in) :: o_mi
    type(t_orb_m), intent(out) :: o_mo
    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: eps, om
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi, &
         drad = Pi/180.d0, mu = TwoPi**2
    type(t_orb_m) :: o_md
    type(t_v3d) :: posi, poso, veli, velo

    o_md%a = o_mi%a
    o_md%e = o_mi%e
    o_md%inc = o_mi%inc
    o_md%node = o_mi%node
    o_md%peri = o_mi%peri
    o_md%m = o_mi%m
    call coord_cart (mu, o_md, posi, veli)
    call ref_ecl (ieqec, posi, poso, eps, om, ierr)
    call ref_ecl (ieqec, veli, velo, eps, om, ierr)
    call osc_el (mu, poso, velo, o_mo)
    call ztopi (o_mo%inc)
    if (o_mo%inc .gt. Pi) o_mo%inc = o_mo%inc - Pi
    call ztopi (o_mo%node)
    call ztopi (o_mo%peri)
    call ztopi (o_mo%m)

    return
  end subroutine ref_ecl_osc

  subroutine ref_ecl_inc_node(ieqec, ii, noi, io, noo, eps, om, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine converts inclination and longitude of node elements back
! and forth between given reference frame and ecliptic plane.
! Uses ref_ecl_osc to do the work.
! Reference frame is given by Omega (om), the node longitude (rotation
! around Z axis of ecliptic) and Epsilon (eps), the inclination of the
! reference frame (rotation around node axis of ecliptic).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) ieqec
!f2py intent(in) ii
!f2py intent(in) noi
!f2py intent(in) eps
!f2py intent(in) om
!f2py intent(out) io
!f2py intent(out) noo
!f2py intent(out) ierr
    implicit none

    integer, intent(in) :: ieqec
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: ii, noi, eps, om
    real (kind=8), intent(out) :: io, noo
    type(t_orb_m) :: o_mi, o_mo

    o_mi%a = 10.d0
    o_mi%e  =0.2d0
    o_mi%peri = 0.d0
    o_mi%m = 0.d0
    o_mi%inc = ii
    o_mi%node = noi
    call ref_ecl_osc (ieqec, o_mi, o_mo, eps, om, ierr)

    return
  end subroutine ref_ecl_inc_node

  subroutine forced_plane(a, ifd, Omfd)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine returns the inclination iforced and node Omforced of the
! forced plane as a function of semimajor-axis in the region of the
! main classical Kuiper belt, according to second order theory with all
! 8 planets.
!
! The fit is valid only between 35 au and 47.74 au (2:1 MMR).
! Returns 0 for a < 35 au, and the invariable plane for a > 47.74 au
!
! Returns degrees.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2018
! Version 2 : March 2019
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     a     : semimajor-axis [au] (R8)
!
! OUTPUT
!     ifd   : Inclination of forced plane [deg] (R8)
!     Omfd  : Node of forced plane [deg] (R8)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in) a
!f2py intent(out) ifd
!f2py intent(out) Omfd
    implicit none

! Calling arguments
    real (kind=8), intent(in) ::  a
    real (kind=8), intent(out) ::  ifd, Omfd

! Internal variables
    real (kind=8) :: ci0(4), ci1(3), ci2(3), ci3(3), co0(3), co1(3), co2(5), &
         co3(5), alpha, beta, gamma, delta, om_lim_low, om_lim_high

    data &
         ci0 /-0.115657583d0,34.8097343d0,7.79198557d-02,-1.06252408d0/, &
         ci1 /-0.392018467d0, 40.5282974d0, 0.137285933/, &
         ci2 /-0.391531527d0, 40.6837921d0, 0.190758109/, &
         ci3 /0.391087621d0, 40.2941360d0, 0.146276891/, &
         co0 /971.35006612732775d0, -50.061061665365997d0, &
              0.74715152013989472d0/, &
         co1 /4694.6567730091410d0, -249.62252879980477d0, &
              3.4214148236506192d0/, &
         co2 /11.7058372d0, -294.139587d0, 30.8862801d0, &
              -1049.41772d0, 15.2022839d0/, &
         co3 /30.4611244d0, -1203.70593d0, 2.47961617d0, &
              -16.6130295d0, 302.805359d0/, &
         alpha /1.d0/, &
         om_lim_low /200.d0/, &
         om_lim_high /20.d0/

    common /om_lim_com/ om_lim_low, om_lim_high

    if (a .lt. 35.d0) then
       ifd = 0.d0
       Omfd = 0.d0
    else if (a .lt. 37.d0) then
       ifd = (ci0(1)/(a - ci0(2)) + ci0(3)*a + ci0(4))
       Omfd = co0(1) + co0(2)*a + co0(3)*a*a
!    else if (a .lt. 39.41d0) then
    else if (a .lt. 39.0d0) then
       ifd = 10.d0**(ci1(1)/(a - ci1(2)) + ci1(3))
       Omfd = co1(1) + co1(2)*a + co1(3)*a*a
    else if (a .lt. 40.47d0) then
!       ifd = 10.d0**(ci2(1)/(a - ci2(2)) + ci2(3))
!       beta = -((co2(1) + co2(3))*a + co2(2) + co2(4))
!       gamma = (co2(1)*a + co2(2))*(co2(3)*a + co2(4)) - co2(5)
!       delta = max(0., beta**2 - 4.d0*alpha*gamma)
!! Here we use the (-b + sqrt(\Delta))/(2a) solution
!       Omfd = (-beta + sqrt(delta))/(2.d0*alpha)
       ifd = 2.477d0
!       Omfd = (200.d0-163.187d0)/(40.47d0-39.d0)*(a-40.47d0)+200.d0
       Omfd = (om_lim_low-163.187d0)/(40.47d0-39.d0)*(a-40.47d0) &
            + om_lim_low
    else if (a .lt. 41.9d0) then
       ifd = 2.477d0
!       Omfd = (61.791d0-20.d0)/(41.9d0-40.47d0)*(a-40.47d0) + 20.d0
       Omfd = (61.791d0-om_lim_high)/(41.9d0-40.47d0)*(a-40.47d0) &
            + om_lim_high
    else if (a .lt. 47.74d0) then
       ifd = 10.d0**(ci3(1)/(a - ci3(2)) + ci3(3))
       beta = -((co3(1) + co3(3))*a + co3(2) + co3(4))
       gamma = (co3(1)*a + co3(2))*(co3(3)*a + co3(4)) - co3(5)
       delta = max(0., beta**2 - 4.d0*alpha*gamma)
! Here we use the (-b - sqrt(\Delta))/(2a) solution
       Omfd = (-beta - sqrt(delta))/(2.d0*alpha)
    else
       ifd = 5713.86d0/3600.d0
       Omfd = 387390.8d0/3600.d0
    end if
    ifd = min(ifd, 40.d0)
 
    return

  end subroutine forced_plane

  subroutine ztopi (var)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This subroutine resets var to be between 0 and 2Pi.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) var
    implicit none

    real (kind=8), intent(inout) :: var
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi

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
  end subroutine ztopi

end module rot
