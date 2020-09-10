module modelutils
  implicit none

  interface
     function func(nparam, param, inc)
       real (kind=8) :: func
       integer, intent(in) :: nparam
       real (kind=8), intent(in) :: param(*), inc
     end function func
  end interface

  type func_holder
     procedure(func), pointer, nopass :: f_ptr => null()
  end type func_holder

contains
  subroutine incdism (seed, nparam, param, incmin, incmax, inc, &
       dist, ierr, func)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a variable according to probability
! density \verb|func| with parameters \verb|param|. Same as previous
! routine, but can remember up to 10 different distributions at a time.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     incmin: Minimum inclination (R8)
!     incmax: Maximum inclination (R8)
!     dist  : index of the selected distribution (I4)
!     func  : probability density function
!
! OUTPUT
!     inc   : Inclination (R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : wrong input data
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) incmin
!f2py intent(in) incmax
!f2py intent(in) dist
!f2py intent(out) inc
!f2py intent(out) ierr
    implicit none

    integer, intent(in) :: nparam, dist
    integer, intent(inout) :: seed
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: param(*), incmin, incmax
    real (kind=8), intent(out) :: inc
    type(func_holder), intent(in) :: func
    integer :: i, ilo, ihi, di
    real (kind=8) :: random
    integer, parameter :: np = 10000, nd = 10
    real (kind=8), save :: proba(0:np,nd), inctab(0:np,nd)
    logical, save :: first(nd)

    data first /.true.,.true.,.true.,.true.,.true., &
         .true.,.true.,.true.,.true.,.true./

    ierr = 0
    di = min(nd, dist)
    if (first(di)) then
       inctab(0,di) = incmin
       proba(0,di) = 0.d0
       do i = 1, np
          inctab(i,di) = incmin + dfloat(i)*(incmax-incmin)/dfloat(np)
          proba(i,di) = func%f_ptr(nparam, param(1:nparam), inctab(i,di)) + proba(i-1,di)
       end do
       do i = 1, np
          proba(i,di) = proba(i,di)/proba(np,di)
       end do
       first(di) = .false.
    end if

    random = ran_3(seed)
    ilo = 0
    ihi = np
1000 continue
    if (ihi - ilo .gt. 1) then
       i = (ihi + ilo)/2
       if (proba(i,di) .lt. random) then
          ilo = i
       else if (proba(i,di) .gt. random) then
          ihi = i
       else
          inc = inctab(i,di)
          return
       end if
       goto 1000
    end if
    inc = inctab(ilo,di) + (inctab(ihi,di) - inctab(ilo,di))* &
         (random - proba(ilo,di))/(proba(ihi,di) - proba(ilo,di))

    return
  end subroutine incdism

  real (kind=8) function offgau (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density as a non-zero centered gaussian.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [Radian] (R8)
!
! OUPUT
!     offgau: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t0, t1, t3

    if (nparam .ne. 2) stop
    t0 = param(1)
    s1 = param(2)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -(t0-angle)**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    offgau = dsin(angle)*fe

    return
  end function offgau

  real (kind=8) function cold_low_a_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [Radian] (R8)
!
! OUPUT
!     cold_low_a_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.
    real (kind=8) :: fe, angle

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. 0.8d0*drad) then
       cold_low_a_inc = angle/(0.8*drad)
    else if (angle .lt. 2.6d0*drad) then
       cold_low_a_inc = 1.d0
    else if (angle .le. 4.3d0*drad) then
       cold_low_a_inc = 0.7d0*(4.3d0*drad-angle)/(1.7d0*drad) + 0.3d0
!    else if (angle .le. 5.d0*drad) then
!       cold_low_a_inc = 0.3d0*(5.d0*drad - angle)/(0.7d0*drad)
    else if (angle .le. 4.7d0*drad) then
       cold_low_a_inc = 0.3d0*(4.7d0*drad - angle)/(0.4d0*drad)
    else
       cold_low_a_inc = 0.d0
    end if
      
    return
  end function cold_low_a_inc

  real (kind=8) function cold_high_a_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [Radian] (R8)
!
! OUPUT
!     cold_high_a_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.
    real (kind=8) :: fe, angle

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. 0.9d0*drad) then
       cold_high_a_inc = angle/(0.9d0*drad)
!    if (angle .lt. 0.8d0*drad) then
!       cold_high_a_inc = angle/(0.8d0*drad)
    else if (angle .le. 2.7d0*drad) then
!    else if (angle .le. 2.6d0*drad) then
       cold_high_a_inc = 1.d0
    else if (angle .le. 4.3d0*drad) then
       cold_high_a_inc = 0.3d0*(4.3d0*drad-angle)/(1.6d0*drad) + 0.7d0
!    else if (angle .le. 9.d0*drad) then
!       cold_high_a_inc = (9.d0*drad - angle)/(6.3d0*drad)
    else if (angle .le. 12.d0*drad) then
       cold_high_a_inc = 0.7d0*(12.d0*drad - angle)/(7.7d0*drad)
    else
       cold_high_a_inc = 0.d0
    end if

    return
  end function cold_high_a_inc

  real (kind=8) function hot_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [Radian] (R8)
!
! OUPUT
!     hot_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.
    real (kind=8) :: fe, angle

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. 10.d0*drad) then
!       hot_inc = sin(angle)
       hot_inc = 0.5d0*angle/(10.d0*drad)
    else if (angle .le. 20.d0*drad) then
       hot_inc = 0.5d0
    else if (angle .le. 50.d0*drad) then
!       hot_inc = 0.5d0*(50.d0*drad - angle)/(20.d0*drad)
       hot_inc = 1.d0 - cos((50.d0*drad-angle)/30.d0*60.d0)
    else
       hot_inc = 0.d0
    end if
 
    return
  end function hot_inc

  real (kind=8) function onecomp (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density of Brown.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2007
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination (R8)
!
! OUPUT
!     onecomp: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t1, t3

    if (nparam .ne. 1) stop
    s1 = param(1)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -angle**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    onecomp = dsin(angle)*fe

    return
  end function onecomp

  real (kind=8) function onecompjmp (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns an inclination drawn according to a single component
! Brown distribution function, modified to include sin(i)**2 instead of sin(i).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2020
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination (R8)
!
! OUPUT
!     onecompjmp: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer, intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t1, t3

    if (nparam .ne. 1) stop
    s1 = param(1)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -angle**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    onecompjmp = dsin(angle)**2*fe

    return
  end function onecompjmp

  real (kind=8) function interp (x, y, val, n)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function linearly interpolates the function y(x) at value x=val.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     x     : Abscissa of the function, sorted in ascending order (n*R8)
!     y     : Values of the function (n*R8)
!     val   : Value of x at which to interpolate (R8)
!     n     : Size of x and y arrays (I4)
!
! OUTPUT
!     interp: Interpolated value (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) n
!f2py intent(in), depend(n) :: x
!f2py intent(in), depend(n) :: y
!f2py intent(in) val
    implicit none

    integer, intent(in) ::  n
    real (kind=8), intent(in) :: x(*), y(*), val
    integer ::  ilo, ihi, i

    if (val .le. x(1)) then
       interp = y(1)
    else if (val .ge. x(n)) then
       interp = y(n)
    else
       ilo = 1
       ihi = n
1000   continue
       if (ihi - ilo .gt. 1) then
          i = (ihi + ilo)/2
          if (x(i) .lt. val) then
             ilo = i
          else if (x(i) .gt. val) then
             ihi = i
          else
             interp = y(i)
             return
          end if
          goto 1000
       end if
       interp = y(ilo) + (y(ihi) - y(ilo))*(val - x(ilo))/(x(ihi) - x(ilo))
    end if

    return
  end function interp

  real (kind=8) function size_dist_one (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to an exponential differential
! distribution specified by "h_params":
!
!   P(h) d_h = A \exp{(h*h_params(3))} d_h
!
! with h_params(1) <= h <= h_params(2)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_one: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer, intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(*)
    real (kind=8) :: h0s10, h1s10, random, slope
!
! H-mag distribution 
!
! Functions have been triple-checked as of 2018-05-04. OK.
    slope = h_params(3)
    h0s10 = 10.d0**(slope*h_params(1))
    h1s10 = 10.d0**(slope*h_params(2))
    random=ran_3(seed)
    size_dist_one = log10( random*(h1s10 - h0s10) + h0s10 ) / slope

    return
  end function size_dist_one

  real (kind=8) function size_dist_two (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to a 2slope exponential
! differential distribution specified by "h_params":
!
!   P(h) d_h = A \exp{(h*h_params(4))} d_h
!
! with h_params(1) <= h <= h_params(2)
!
!   P(h) d_h = B \exp{(h*h_params(5))} d_h
!
! with h_params(2) <= h <= h_params(3)
!
! continuous at H1 = h_params(2)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_two: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer, intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(5)
    real (kind=8) :: h0s1, h1s1, h1s2, h2s2, h1s12, random, sl1, sl2, xi1
!
! H-mag distribution
!
! Functions have been triple-checked as of 2018-05-04. OK.
    sl1 = h_params(4)
    sl2 = h_params(5)
    h0s1 = 10.d0**(sl1*h_params(1))
    h1s1 = 10.d0**(sl1*h_params(2))
    h1s2 = 10.d0**(sl2*h_params(2))
    h1s12 = h1s1/h1s2
    h2s2 = 10.d0**(sl2*h_params(3))
    xi1 = sl2*(h1s1-h0s1)/(sl1*h1s12*(h2s2-h1s2)+sl2*(h1s1-h0s1))
    random=ran_3(seed)
    if (random .le. xi1) then
       size_dist_two = log10(random/xi1*(h1s1 - h0s1) + h0s1)/sl1
    else
       size_dist_two = log10((random-xi1)/(1.d0-xi1)*(h2s2 - h1s2) &
            + h1s2)/sl2
    end if

    return
  end function size_dist_two

  real (kind=8) function size_dist_two_straight (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to a 2-slope exponential
! cumulative distribution specified by "h_params" (2 straight lines in
! semilog for cumulative):
!
!   P(<=h) = A 10^{(h*h_params(4))}
!
! with h_params(1) <= h <= h_params(2)
!
!   P(<=h) = B 10^{(h*h_params(5))}
!
! with h_params(2) < h <= h_params(3)
!
! continuous at H1 = h_params(2)
!   (A = B 10^{(h_params(5)-h_params(4))*h_params(2)})
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : May 2018
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_two_straight: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer, intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(5)
    real (kind=8) :: h0s1, h1s1, h1s2, h2s2, h12s2, random, sl1, sl2
!
! H-mag distribution
!
! Functions have been triple-checked as of 2018-05-04. OK.
!
! In order to not have the ramp-up found in {\it size_dist_two}, I use
! an exponential with no lower limit on H. What sets the lower limit of
! the faint end is the selection of the correct slope at H_break (see
! test below). The test is based on the number of objects in each
! portion of the distribution.
!
! Since we want N(<H) = A 10^{\alpha_1 H} for H < H_b
! and           N(<H) = B 10^{\alpha_2 H} for H >= H_b
! and N(<H) continuous at H_b, we have N(<H_b) = A 10^{\alpha_1 H_b} =
! B 10^{\alpha_2 H_b}.
!
! The total number of objects is N(<H_max) = B 10^{\alpha_2 H_max}, and
! the number of objects bigger than H_b is N(<H_b) = B 10^{\alpha_2 H_b}. 
! So the fraction of objects bigger than H_b is 10^{\alpha_2 H_b} /
! 10^{\alpha_2 H_max} = h1s2 / h2s2 = h12s2.
!
! Therefore, if 'random' < h12s2, then we have a big object, and
! otherwise a small on. 'random' needs to be rescaled to go to 1.
    sl1 = h_params(4)
    sl2 = h_params(5)
    h0s1 = 10.d0**(sl1*h_params(1))
    h1s1 = 10.d0**(sl1*h_params(2))
    h1s2 = 10.d0**(sl2*h_params(2))
    h2s2 = 10.d0**(sl2*h_params(3))
    h12s2 = h1s2/h2s2
    random=ran_3(seed)
    if (random .le. h12s2) then
       size_dist_two_straight = log10(random*h1s1/h12s2)/sl1
    else
       size_dist_two_straight = log10(random*h2s2)/sl2
    end if

    return
  end function size_dist_two_straight

  real (kind=8) function size_dist_n_straight (seed, h_params, nn)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to an n-slope exponential
! cumulative distribution specified by "h_params" (n straight lines in
! semilog for cumulative):
!
!   P(<=h) = A_k 10^{(h*\alpha_k)}
!
! with H_{k-1} < h <= H_k,
!
! with
!
! \alpha_k = h_params(k+n)
! H_k = h_params(k)
!
! for k in [1; n]. h_params(0) is not present and implicit at -\infty.
! This corresponds to dropping the lower limit as we want straight lines.
!
! The function is continuous at H_k = h_params(k) ofr all k's:
!   A_k 10^{H_k*\alpha_k} = A_{k+1} 10^{H_k*\alpha_{k+1}}
!
! for k in [1; n-1].
!
! To avoid using allocatale arrays and dynamical allocation, I restrict
! the number of slopes to be <= 10.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2019
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!     nn    : number of different slopes (I4)
!
! OUTPUT
!     size_dist_n_straight: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
!f2py intent(in) nn
    implicit none

    integer, intent(in) :: nn
    integer, intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(*)
    integer :: n, k
    real (kind=8) :: random, x(10), f(10)
    logical, save :: first

    data first /.true./

    if (nn .gt. 10) then
       if (first) then
          first = .false.
          print *, &
               'WARNING: number of breaks greater than allowed maximum'// &
               ' 10, using only the 10 first values.'
       end if
       n = 10
    else
       n = nn
    end if
!
! H-mag distribution
!
! In order to not have the ramp-up found in {\it size_dist_two}, I use
! an exponential with no lower limit on H. What sets the lower limit of
! the faint end is the selection of the correct slope at H_k (see
! test below). The test is based on the number of objects in each
! portion of the distribution.
!
! The total number of objects is N(<H_max) = N(<H_n) = A_n 10^{H_n \alpha_n},
! and the number of objects bigger than the next break H_{n-1} is N(
! <H_{n-1}) = A_n 10^{H_{n-1} \alpha_n}. More generaly, the number of
! objects bigger than H_k is
!
! N_k = A_k 10^{H_k \alpha_k} = A_{k+1} 10^{H_k \alpha_{k+1}}
!
! Let's define the fraction of objects bigger than H_k in the objects
! bigger than H_{k+1} as:
!
! x_k = N_{k-1} / N_k = 10^{(H_{k-1}-H_k) \alpha_k}
!
! for k in [2; n], and x_1 = 0.
!
! The fraction of object bigger than B_k compared to the total
! number of objects is
!
! f_k = N_k / N_n = (N_k / N_{k+1}) (N_{k+1} / N_{k+2}) ... (N_{n-1} / N_n)
!                 = x_{k+1} x_{k+2} ... x_n
!
! We define x_k as described above, then we set
!
! f_n = 1.
! f_{k-1} = x_k f_k, for k in [2; n]
!
! Therefore, if (f_{k-1} < random <= f_k), then we have an object in
! range ]H_{k-1}; H_k]. 'random' needs to be rescaled to go to 1.
    x(1) = 0.d0
    do k = 2, n
       x(k) = 10.d0**((h_params(k-1)-h_params(k))*h_params(k+n))
    end do
    f(n) = 1.d0
    do k = n, 2, -1
       f(k-1) = f(k)*x(k)
    end do
    random=ran_3(seed)
    do k = 1, n
       if (random .le. f(k)) then
          size_dist_n_straight = &
               log10(random/f(k)*10.d0**(h_params(k)*h_params(k+n))) &
               /h_params(k+n)
          return
       end if
    end do

  end function size_dist_n_straight

  real (kind=8) function qhot (np, p, q)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This function returns a 'q' drawn according to the probability of having
! perihelion distance "q"
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2020
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     np    : Number of parameters describing the distribution function (I4)
!     p     : Array of parameters (np*R8)
!     q     : Perihelion distance (R8)
!
! OUTPUT
!     qhot  : Probability of q
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) q
    implicit none

! Calling arguments
    integer, intent(in) :: np
    real (kind=8), intent(in) :: p(*), q

    qhot = 1./((1.+exp((p(3)-q)/p(4)))*(1.+exp((q-p(1))/p(2))))

    return
  end function qhot

  real (kind=8) function ran_3(idum)
!f2py intent(in,out) idum
    INTEGER, intent(inout) :: idum
    INTEGER, parameter :: MBIG=1000000000, MSEED=161803398, MZ=0
    REAL (KIND=8), parameter :: FAC=1.d0/MBIG
    INTEGER :: i,ii,k,mj,mk
    INTEGER, save :: iff,inext,inextp,ma(55)
    data iff /0/
    if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=abs(MSEED-abs(idum))
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
       end do
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
       end do
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran_3=mj*FAC
    return
  end function ran_3

end module modelutils
