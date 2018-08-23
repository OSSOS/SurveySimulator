      subroutine psalun (i, x)
c********************************************************************
c
c This subroutine computes a "random" number using the recurrence
c formala:
c i(k+1) = i(k) * 367379597 + 1 mod(2**31).
c
c x is a double precision real "random" number between 0 (inclusive)
c and 1 (exclusive). x = i/2.**31
c
c Upon first call to psalun, one must give an initial value to i
c (a seed) and then never change i (updated by the subroutine).
c
c********************************************************************
      integer*4
     1  i, k, mask

      real*8
     1  x

      parameter
     1  (k = 367379597)

      data
     1  mask /x'7fffffff'/

      i = i*k + 1
      i = iand(i,mask)
      x = dfloat(i)/2147483648.d0

      return
      end

      subroutine incdis (seed, nparam, param, incmin, incmax, inc,
     $  ierr, func)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly an inclination according to probability
c density \verb|func| with parameters \verb|param|.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     incmin: Minimum inclination (R8)
c     incmax: Maximum inclination (R8)
c     func  : probability density function
c
c OUTPUT
c     inc   : Inclination (R8)
c     ierr  : Error code
c                0 : nominal run
c               10 : wrong input data
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) incmin
Cf2py intent(in) incmax
Cf2py intent(out) inc
Cf2py intent(out) ierr

      implicit none

      integer
     $  np

      parameter
     $  (np = 10000)

      integer
     $  ierr, seed, nparam, i, ilo, ihi

      real*8
     $  param(*), inc, proba(0:np), inctab(0:np), random, func,
     $  incmin, incmax, ran3

      logical
     $  first

      external
     $  func, ran3

      save proba, inctab, first

      data first /.true./

      ierr = 0
      if (first) then
         inctab(0) = incmin
         proba(0) = 0.d0
         do i = 1, np
            inctab(i) = incmin + dfloat(i)*(incmax-incmin)/dfloat(np)
            proba(i) = func(nparam, param, inctab(i)) + proba(i-1)
         end do
         do i = 1, np
            proba(i) = proba(i)/proba(np)
         end do
         first = .false.
      end if

      random = ran3(seed)
      ilo = 0
      ihi = np
 1000 continue
      if (ihi - ilo .gt. 1) then
         i = (ihi + ilo)/2
         if (proba(i) .lt. random) then
            ilo = i
         else if (proba(i) .gt. random) then
            ihi = i
         else
            inc = inctab(i)
            return
         end if
         goto 1000
      end if
      inc = inctab(ilo) + (inctab(ihi) - inctab(ilo))*
     $  (random - proba(ilo))/(proba(ihi) - proba(ilo))

      return
      end

      subroutine incdism (seed, nparam, param, incmin, incmax, inc,
     $  dist, ierr, func)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly a variable according to probability
c density \verb|func| with parameters \verb|param|. Same as previous
c routine, but can remember up to 10 different distributions at a time.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     incmin: Minimum inclination (R8)
c     incmax: Maximum inclination (R8)
c     dist  : index of the selected distribution (I4)
c     func  : probability density function
c
c OUTPUT
c     inc   : Inclination (R8)
c     ierr  : Error code
c                0 : nominal run
c               10 : wrong input data
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) incmin
Cf2py intent(in) incmax
Cf2py intent(in) dist
Cf2py intent(out) inc
Cf2py intent(out) ierr

      implicit none

      integer
     $  np, nd

      parameter
     $  (np = 10000, nd = 10)

      integer
     $  ierr, seed, nparam, i, ilo, ihi, dist, di

      real*8
     $  param(*), inc, proba(0:np,nd), inctab(0:np,nd), random, func,
     $  incmin, incmax, ran3

      logical
     $  first(nd)

      external
     $  func, ran3

      save proba, inctab, first

      data first /.true.,.true.,.true.,.true.,.true.,
     $            .true.,.true.,.true.,.true.,.true./

      ierr = 0
      di = min(nd, dist)
      if (first(di)) then
         inctab(0,di) = incmin
         proba(0,di) = 0.d0
         do i = 1, np
            inctab(i,di) = incmin + dfloat(i)*(incmax-incmin)/dfloat(np)
            proba(i,di) = func(nparam, param, inctab(i,di))
     $        + proba(i-1,di)
         end do
         do i = 1, np
            proba(i,di) = proba(i,di)/proba(np,di)
         end do
         first(di) = .false.
      end if

      random = ran3(seed)
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
      inc = inctab(ilo,di) + (inctab(ihi,di) - inctab(ilo,di))*
     $  (random - proba(ilo,di))/(proba(ihi,di) - proba(ilo,di))

      return
      end

      real*8 function offgau (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density as a non-zero centered gaussian.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : August 2014
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c OUPUT
c     offgau: Value of the probability (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc

      real*8
     $  Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi)

      real*8
     $  fe, s1, angle, t0, t1, t3

      if (nparam .ne. 2) stop
      t0 = param(1)
      s1 = param(2)
      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      t1 = 2.*s1**2
      t3 = -(t0-angle)**2
      fe = exp(t3/t1)
      offgau = dsin(inc)*fe

      return
      end

      real*8 function onecomp (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density of Brown.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2007
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c OUPUT
c     onecomp: Value of the probability (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc

      real*8
     $  Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi)

      real*8
     $  fe, s1, angle, t1, t3

      if (nparam .ne. 1) stop
      s1 = param(1)
      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      t1 = 2.*s1**2
      t3 = -angle**2
      fe = exp(t3/t1)
      onecomp = dsin(inc)*fe

      return
      end

      real*8 function onecompjmp (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density of Brown, modified to include sin(i)**2 instead of sin(i).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2007
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c OUPUT
c     onecompjmp: Value of the probability (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc

      real*8
     $  Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi)

      real*8
     $  fe, s1, angle, t1, t3

      if (nparam .ne. 1) stop
      s1 = param(1)
      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      t1 = 2.*s1**2
      t3 = -angle**2
      fe = exp(t3/t1)
      onecompjmp = dsin(inc)**2*fe

      return
      end

      real*8 function brown (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density of Brown.
c This is just a wrapper to brownjmp.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc, brownjmp

      external
     $  brownjmp

      brown = brownjmp (nparam, param, inc)

      return
      end

      real*8 function brownjmp (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density of Brown.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc

      real*8
     $  fe, a, s1, s2

      if (nparam .ne. 3) stop
      a = param(1)
      s1 = param(2)
      s2 = param(3)
      call incecl1 (inc, fe)
      brownjmp = dsin(inc)*fe

      return
      end

      subroutine incecl1 (inc, fe)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the ecliptic inclination distribution for a
c given inclination (given in radian):
c
c   f_e(i) = a exp{-i^2/(2 s_1^2)} + (1 - a) exp{-i^2/(2 s_2^2)}
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : April 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     inc   : inclination requested (radian) (R8)
c
c OUTPUT
c     fe    : ecliptic inclination distribution (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) inc
Cf2py intent(out) fe

      implicit none

      real*8
     $  Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi)

      real*8
     $  inc, fe, a, s1, s2, angle, t1, t2, t3

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      t1 = 2.*s1**2
      t2 = 2.*s2**2
      t3 = -angle**2
      fe = a*exp(t3/t1) + (1. - a)*exp(t3/t2)

      return
      end

      real*8 function interp (x, y, val, n)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function linearly interpolates the function y(x) at value x=val.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     x     : Abscissa of the function, sorted in ascending order (n*R8)
c     y     : Values of the function (n*R8)
c     val   : Value of x at which to interpolate (R8)
c     n     : Size of x and y arrays (I4)
c
c OUTPUT
c     interp: Interpolated value (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) n
Cf2py intent(in), depend(n) :: x
Cf2py intent(in), depend(n) :: y
Cf2py intent(in) val

      implicit none

      integer*4
     $  n, ilo, ihi, i

      real*8
     $  x(*), y(*), val

      if (val .le. x(1)) then
         interp = y(1)
      else if (val .ge. x(n)) then
         interp = y(n)
      else
         ilo = 1
         ihi = n
 1000    continue
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
         interp = y(ilo) + (y(ihi) - y(ilo))*(val - x(ilo))
     $     /(x(ihi) - x(ilo))
      end if

      return
      end

      real*8 function size_dist_one (seed, h_params)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function draws a number according to an exponential differential
c distribution specified by "h_params":
c
c   P(h) d_h = A \exp{(h*h_params(3))} d_h
c
c with h_params(1) <= h <= h_params(2)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : April 2004
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : seed for the random number generator (I4)
c     h_params: parameters for the distribution (3*R8)
c
c OUTPUT
c     size_dist_one: Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in,out) seed
Cf2py intent(in) h_param
c
      implicit none

      integer*4 seed
      real*8 h_params(*), h0s10, h1s10, random, slope, ran3

      external
     $  ran3
c
c H-mag distribution 
c
c Functions have been triple-checked as of 2018-05-04. OK.
      slope = h_params(3)
      h0s10 = 10.d0**(slope*h_params(1))
      h1s10 = 10.d0**(slope*h_params(2))
      random=ran3(seed)
      size_dist_one = log10( random*(h1s10 - h0s10) + h0s10 ) / slope

      return

      end

      real*8 function size_dist_two (seed, h_params)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function draws a number according to a 2slope exponential
c differential distribution specified by "h_params":
c
c   P(h) d_h = A \exp{(h*h_params(4))} d_h
c
c with h_params(1) <= h <= h_params(2)
c
c   P(h) d_h = B \exp{(h*h_params(5))} d_h
c
c with h_params(2) <= h <= h_params(3)
c
c continuous at H1 = h_params(2)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : April 2014
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : seed for the random number generator (I4)
c     h_params: parameters for the distribution (3*R8)
c
c OUTPUT
c     size_dist_two: Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in,out) seed
Cf2py intent(in) h_param
c
      implicit none

      integer*4
     $  seed

      real*8
     $  h_params(5), h0s1, h1s1, h1s2, h2s2, h1s12, random, sl1, sl2,
     $  xi1, ran3

      external
     $  ran3
c
c H-mag distribution
c
c Functions have been triple-checked as of 2018-05-04. OK.
      sl1 = h_params(4)
      sl2 = h_params(5)
      h0s1 = 10.d0**(sl1*h_params(1))
      h1s1 = 10.d0**(sl1*h_params(2))
      h1s2 = 10.d0**(sl2*h_params(2))
      h1s12 = h1s1/h1s2
      h2s2 = 10.d0**(sl2*h_params(3))
      xi1 = sl2*(h1s1-h0s1)/(sl1*h1s12*(h2s2-h1s2)+sl2*(h1s1-h0s1))
      random=ran3(seed)
      if (random .le. xi1) then
         size_dist_two = log10(random/xi1*(h1s1 - h0s1) + h0s1)/sl1
      else
         size_dist_two = log10((random-xi1)/(1.d0-xi1)*(h2s2 - h1s2)
     $     + h1s2)/sl2
      end if

      return

      end

      real*8 function size_dist_two_straight (seed, h_params)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function draws a number according to a 2-slope exponential
c cumulative distribution specified by "h_params" (2 straight lines in
c semilog for cumulative):
c
c   P(<=h) = A 10^{(h*h_params(4))}
c
c with h_params(1) <= h <= h_params(2)
c
c   P(<=h) = B 10^{(h*h_params(5))}
c
c with h_params(2) < h <= h_params(3)
c
c continuous at H1 = h_params(2)
c   (A = B 10^{(h_params(5)-h_params(4))*h_params(2)})
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2018
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : seed for the random number generator (I4)
c     h_params: parameters for the distribution (3*R8)
c
c OUTPUT
c     size_dist_two_straight: Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in,out) seed
Cf2py intent(in) h_param
c
      implicit none

      integer*4
     $  seed

      real*8
     $  h_params(5), h0s1, h1s1, h1s2, h2s2, h12s2, random, sl1, sl2,
     $  ran3

      external
     $  ran3
c
c H-mag distribution
c
c Functions have been triple-checked as of 2018-05-04. OK.
      sl1 = h_params(4)
      sl2 = h_params(5)
      h0s1 = 10.d0**(sl1*h_params(1))
      h1s1 = 10.d0**(sl1*h_params(2))
      h1s2 = 10.d0**(sl2*h_params(2))
      h2s2 = 10.d0**(sl2*h_params(3))
      h12s2 = h1s2/h2s2
      random=ran3(seed)
      if (random .le. h12s2) then
         size_dist_two_straight = log10(random*h1s1/h12s2)/sl1
      else
         size_dist_two_straight = log10(random*h2s2)/sl2
      end if

      return

      end

      real*8 function qhot (np, p, q)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c This function returns the probability of having perihelion distance "q"
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2013
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT
c     np    : Number of parameters describing the distribution function (I4)
c     p     : Array of parameters (np*R8)
c     q     : Perihelion distance (R8)
c
c OUTPUT
c     qhot  : Probability of q
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in) np
Cf2py intent(in) p
Cf2py intent(in) q
c
      implicit none

c Calling arguments
      integer np
      real*8 p(*), q

      qhot = 1./((1.+exp((p(3)-q)/p(4)))*(1.+exp((q-p(1))/p(2))))

      return
      end
