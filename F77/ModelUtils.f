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
     1  mask /z'7fffffff'/

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
     $  ierr, seed, nparam, i

      real*8
     $  param(*), inc, proba(0:np), inctab(0:np), random, func,
     $  incmin, incmax, ran_3, interp

      logical
     $  first

      external
     $  func, ran_3, interp

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

      random = ran_3(seed)
      inc = interp(proba, inctab, random, np+1)

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
     $  ierr, seed, nparam, i, dist, di

      real*8
     $  param(*), inc, proba(0:np,nd), inctab(0:np,nd), random, func,
     $  incmin, incmax, ran_3, interp

      logical
     $  first(nd)

      external
     $  func, ran_3, interp

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

      random = ran_3(seed)
      inc = interp(proba(0,di), inctab(0,di), random, np+1)

      return
      end

      real*8 function Variably_tapered(h, params)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the number of objects brighter or equal to H
c following an exponentially tapered exponential with parameters in
c params.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2021
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     params: parameters for the distribution (4*R8)
c
c OUTPUT
c     Variably_tapered : Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) h
Cf2py intent(in) params

      implicit none

      real*8 params(4), h

      Variably_tapered = 10.d0**(params(3)*3.d0*(h-params(1))/5.d0)
     $     *exp(-10.d0**(-params(4)*3.d0*(h-params(2))/5.d0))

      return
      end

      real*8 function H_dist_cold(seed, h_max)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly a number according to the cold belt H_r
c distribution, represented by an exponentially tapered exponential,
c with parameters frmo the exponential cutoff paper.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2021
c Version 2 : December 2021- updated parameter values.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     h_max : Maximum value of H (R8)
c
c OUTPUT
c     H_dist_cold : Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed

      implicit none

      integer
     $  np

      parameter
     $  (np = 16384)

      integer
     $  seed, nparam, i

      real*8
     $  params(4), proba(0:np), htab(0:np), random,
     $  Variably_tapered, h_min, h_max, ran_3, interp

      logical
     $  first

      external
     $  Variably_tapered, ran_3, interp

      save proba, htab, first

      data first /.true./
      data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
c      data params /-2.466d0, 7.895d0, 0.667d0, 0.438d0/
      data h_min /4.6d0/

      if (first) then
         htab(0) = h_min
         proba(0) = 1.d-10
         do i = 1, np
            htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
            proba(i) = Variably_tapered(htab(i), params)
         end do
         do i = 0, np
            proba(i) = proba(i)/proba(np)
         end do
         first = .false.
      end if

      random = ran_3(seed)
      H_dist_cold = interp(proba, htab, random, np+1)

      return
      end

      real*8 function H_dist_cold_2(seed, nparam, hparam)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly a number according to the cold belt H_r
c distribution, represented by an exponentially tapered exponential,
c with parameters frmo the exponential cutoff paper.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2021
c Version 2 : December 2021 - updated parameter values.
c Version 3 : January 2022 - forcing slope at small sizes.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     nparam: Number of parameters (I4)
c     hparam: Parameters for asymptotic slope(s) (n*R8)
c             hparam(1): start of asymptote
c             hparam(2): contrast at start of asymptote
c             hparam(3): slope of asymptote
c             hparam(4): end of asymptote
c
c OUTPUT
c     H_dist_cold_2 : Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed
      implicit none

      integer
     $  np

      parameter
     $  (np = 16384)

      integer
     $  seed, nparam, i

      real*8
     $  params(4), proba(0:np), htab(0:np), random, hparam(*),
     $  Variably_tapered, h_min, h_max, ran_3, interp, n, c

      logical
     $  first

      external
     $  Variably_tapered, ran_3, interp

      save proba, htab, first

      data first /.true./
      data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
c      data params /-2.466d0, 7.895d0, 0.667d0, 0.438d0/
      data h_min /4.6d0/

      if (first) then
         h_max = hparam(nparam)
         htab(0) = h_min
         proba(0) = 1.d-10
         n = Variably_tapered(hparam(1), params)
         c = hparam(2)
         do i = 1, np
            htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
            if (htab(i) .lt. hparam(1)) then
               proba(i) = Variably_tapered(htab(i), params)
            else
               proba(i) = n
     $           + n*c*(10.0d0**(hparam(3)*(htab(i)-hparam(1))) -
     $           1.0d0)
            end if
         end do
         do i = 0, np
            proba(i) = proba(i)/proba(np)
         end do
         first = .false.
      end if

      random = ran_3(seed)
      H_dist_cold_2 = interp(proba, htab, random, np+1)

      return
      end

      real*8 function H_dist_hot(seed, h_max)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly a number according to the hot belt H_r
c distribution, represented by an exponentially tapered exponential,
c with parameters I've fitted on the OSSOS hot belt data, in range 6
c -8.3. This is the actual maximu likelihood, not the MCMC value.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2021
c Version 2 : December 2021- updated parameter values.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     h_max : Maximum value of H (R8)
c
c OUTPUT
c     H_dist_hot : Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed

      implicit none

      integer
     $  np

      parameter
     $  (np = 16384)

      integer
     $  seed, nparam, i

      real*8
     $  params(4), proba(0:np), htab(0:np), random,
     $  Variably_tapered, h_min, h_max, ran_3, interp, n1, n2, n3,
     $  sl1, sl2, sl3, c1, cb, c, h1, h2, h3

      logical
     $  first

      external
     $  Variably_tapered, ran_3, interp

      save proba, htab, first

      data first /.true./
      data params /-2.465d0, 7.114d0, 0.666d0, 0.875d0/
      data h_min /-1.d0/
      data sl1 /0.13d0/, sl2 /0.5d0/, sl3 /0.4d0/, c /1.d0/
      data h1 /3.2d0/, h2 /6.d0/, h3 /8.5d0/, n1 /3.d0/

      if (first) then
         htab(0) = h_min
         proba(0) = 1.d-10
         n2 = Variably_tapered(h2, params)
         n3 = Variably_tapered(h3, params)
c The normalisation is done with the exponentially tapered exponential,
c as fitted on the OSSOS hot component. This determines the normalisation
c of the exponential between H = h2 and H = h1
c N(<H) = n2*10**(sl2*(H-h2))
c Then, there is an excess divot at h1. See
c [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
c for the appropriate formula.
         cb = n1*sl1*log(10.d0)
         c1 = (n2-n1)*sl2/(n1*sl1*(10.d0**(sl2*(h2-h1))-1.d0))
         do i = 1, np
            htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
            if (htab(i) .lt. h1) then
               proba(i) = n1*10.d0**(sl1*(htab(i)-h1))
            else if (htab(i) .lt. h2) then
               proba(i) = n1 + c1*cb
     $           *(10.d0**(sl2*(htab(i)-h1))-1.d0)/(sl2*log(10.d0))
            else if (htab(i) .lt. h3) then
               proba(i) = Variably_tapered(htab(i), params)
            else
               proba(i) = n3 + n3*c*(10.d0**(sl3*(htab(i)-h3)) - 1.d0)
            end if
         end do
         do i = 0, np
            proba(i) = proba(i)/proba(np)
         end do
         first = .false.
      end if

      random = ran_3(seed)
      H_dist_hot = interp(proba, htab, random, np+1)

      return
      end

      real*8 function H_dist_hot_2(seed, nparam, hparam)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine draws randomly a number according to the hot belt H_r
c distribution, represented by an exponentially tapered exponential,
c with parameters I've fitted on the OSSOS hot belt data, in range 6
c -8.3. This is the actual maximu likelihood, not the MCMC value.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2021
c Version 2 : December 2021- updated parameter values.
c Version 3 : January 2022 - forcing slope at small sizes.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : Random number generator seed (I4)
c     nparam: Number of parameters (I4)
c     hparam: Parameters for asymptotic slope(s) (n*R8)
c             hparam(1): start of asymptote
c             hparam(2): contrast at start of asymptote
c             hparam(3): slope of asymptote
c             hparam(4): end of asymptote
c
c OUTPUT
c     H_dist_hot_2 : Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in,out) seed

      implicit none

      integer
     $  np

      parameter
     $  (np = 16384)

      integer
     $  seed, nparam, i

      real*8
     $  params(4), proba(0:np), htab(0:np), random, hparam(*),
     $  Variably_tapered, h_min, h_max, ran_3, interp, n1, n2, n3,
     $  sl1, sl2, sl3, c1, cb, c, h1, h2, h3

      logical
     $  first

      external
     $  Variably_tapered, ran_3, interp

      save proba, htab, first

      data first /.true./
      data params /-2.465d0, 7.114d0, 0.666d0, 0.875d0/
      data h_min /-1.d0/
      data sl1 /0.13d0/, sl2 /0.5d0/
      data h1 /3.2d0/, h2 /6.d0/, n1 /3.d0/

      if (first) then
         h_max = hparam(nparam)
         htab(0) = h_min
         proba(0) = 1.d-10
         h3 = hparam(1)
         n2 = Variably_tapered(h2, params)
         n3 = Variably_tapered(h3, params)
         c = hparam(2)
         sl3 = hparam(3)
c The normalisation is done with the exponentially tapered exponential,
c as fitted on the OSSOS hot component. This determines the normalisation
c of the exponential between H = h2 and H = h1
c N(<H) = n2*10**(sl2*(H-h2))
c Then, there is an excess divot at h1. See
c [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
c for the appropriate formula.
         cb = n1*sl1*log(10.d0)
         c1 = (n2-n1)*sl2/(n1*sl1*(10.d0**(sl2*(h2-h1))-1.d0))
         do i = 1, np
            htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
            if (htab(i) .lt. h1) then
               proba(i) = n1*10.d0**(sl1*(htab(i)-h1))
            else if (htab(i) .lt. h2) then
               proba(i) = n1 + c1*cb
     $           *(10.d0**(sl2*(htab(i)-h1))-1.d0)/(sl2*log(10.d0))
            else if (htab(i) .lt. h3) then
               proba(i) = Variably_tapered(htab(i), params)
            else
               proba(i) = n3 + n3*c*(10.d0**(sl3*(htab(i)-h3)) - 1.d0)
            end if
         end do
         do i = 0, np
            proba(i) = proba(i)/proba(np)
         end do
         first = .false.
      end if

      random = ran_3(seed)
      H_dist_hot_2 = interp(proba, htab, random, np+1)

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
c     inc   : Inclination [rad] (R8)
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

      real*8 function cold_low_a_inc (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density for cold objects with a < 44.4
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : December 2019
c Version 2 : December 2019, 20th
c Version 3 : December 2019, 28th
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination [rad] (R8)
c
c OUPUT
c     cold_low_a_inc: Value of the probability (R8)
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
     $  Pi, TwoPi, drad

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, drad = Pi/180.)
      real*8 a1, a2, a3, a4
      parameter (a1 = 0.8d0, a2 = 2.6d0, a3 = 4.3d0, a4 = 4.7d0)

      real*8
     $  angle

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      if (angle .lt. a1*drad) then
         cold_low_a_inc = angle/(a1*drad)
      else if (angle .lt. a2*drad) then
         cold_low_a_inc = 1.d0
      else if (angle .le. a3*drad) then
         cold_low_a_inc = 0.7d0*(a3*drad-angle)/((a3-a2)*drad) + 0.3d0
      else if (angle .le. a4*drad) then
         cold_low_a_inc = 0.3d0*(a4*drad - angle)/((a4-a3)*drad)
      else
         cold_low_a_inc = 0.d0
      end if

      return
      end

      real*8 function cold_low_a_inc_2 (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density for cold objects with a < 44.4
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : December 2019
c Version 2 : December 2019, 20th
c Version 3 : December 2019, 28th
c Version 4 : December 2021, 6th
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination [rad] (R8)
c
c OUPUT
c     cold_low_a_inc_2: Value of the probability (R8)
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
     $  Pi, TwoPi, drad

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, drad = Pi/180.)
      real*8 a1, a2, a3
      parameter (a1 = 0.3d0, a2 = 4.05d0, a3 = 4.5d0)

      real*8
     $  angle

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      if (angle .lt. a1*drad) then
         cold_low_a_inc_2 = angle/(a1*drad)
      else if (angle .lt. a2*drad) then
         cold_low_a_inc_2 = 1.d0
      else if (angle .le. a3*drad) then
         cold_low_a_inc_2 = (a3*drad - angle)/((a3-a2)*drad)
      else
         cold_low_a_inc_2 = 0.d0
      end if

      return
      end

      real*8 function kernel_inc (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density for cold objects in the kernel region.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : December 2021, 7th
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination [rad] (R8)
c
c OUPUT
c     onecomp: Value of the probability (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer nparam

      real*8 param(*), inc

      real*8 Pi, TwoPi, drad

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, drad = Pi/180.)
      real*8 a1, a2, delta
      parameter (a1 = 0.9d0, a2 = 2.0d0, delta = 0.8d0)

      real*8 angle, fkrla, fkla, fli, sum_la, da, lli, lhi
      real*8 cold_low_a_inc_2
      integer i, nstep
      logical first

      common /com_kernel/ fkrla, fkla, fli

      external cold_low_a_inc_2
      data first /.true./, nstep /5000/

      save

      if (first) then
         sum_la = 0.d0
         da = 5.d0*drad/dble(nstep)
         do i = 1, nstep
            angle = i*da
            sum_la = sum_la + cold_low_a_inc_2(nparam, param, angle)
         end do
         sum_la = sum_la*da*fkla/(fkrla-fkla)
         lli = fli*sum_la/(delta*drad)
         lhi = (1.d0-fli)*sum_la/(delta*drad)
      end if

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      kernel_inc = cold_low_a_inc_2(nparam, param, angle)
      if ((angle .ge. a1*drad) .and. (angle .le. (a1+delta)*drad)) then
         kernel_inc = kernel_inc + lli
      else if ((angle .ge. a2*drad).and.(angle .le. (a2+delta)*drad))
     $     then
         kernel_inc = kernel_inc + lhi
      end if

      return
      end

      real*8 function cold_high_a_inc (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density for cold objects with a < 44.4
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : December 2019
c Version 2 : December 2019, 20th
c Version 3 : December 2019, 28th
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination [rad] (R8)
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
     $  Pi, TwoPi, drad

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, drad = Pi/180.)
      real*8 a1, a2, a3, a4
      parameter (a1 = 1.2d0, a2 = 4.0d0, a3 = 6.0d0, a4 = 12.0d0)

      real*8
     $  angle

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      if (angle .lt. a1*drad) then
         cold_high_a_inc = angle/(a1*drad)
      else if (angle .le. a2*drad) then
         cold_high_a_inc = 1.d0
      else if (angle .le. a3*drad) then
         cold_high_a_inc = 0.3d0*(a3*drad-angle)/((a3-a2)*drad) + 0.7d0
      else if (angle .le. a4*drad) then
         cold_high_a_inc = 0.7d0*(a4*drad - angle)/((a4-a3)*drad)
      else
         cold_high_a_inc = 0.d0
      end if

      return
      end

      real*8 function hot_inc (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized inclination "probability"
c density for cold objects with a < 44.4
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : December 2019
c Version 2 : December 2019, 20th
c Version 3 : December 2019, 28th
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination [rad] (R8)
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
     $  Pi, TwoPi, drad

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi, drad = Pi/180.)
      real*8 a0, a1, a2, a3
      parameter (a0 = 5.0d0, a1 = 7.0d0, a2 = 26.0d0, a3 = 46.0d0)

      real*8
     $  angle

      angle = mod(inc, TwoPi)
      if (angle .gt. Pi) angle = angle - TwoPi
      if (angle .lt. a0*drad) then
         hot_inc = 0.d0
      else if (angle .lt. a1*drad) then
         hot_inc = 0.5d0*(angle-a0*drad)/((a1-a0)*drad)
      else if (angle .le. a2*drad) then
         hot_inc = 0.5d0
      else if (angle .le. a3*drad) then
         hot_inc = 1.d0 - cos((a3*drad-angle)/(a3-a2)*60.d0)
      else
         hot_inc = 0.d0
      end if

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
      real*8 h_params(*), h0s10, h1s10, random, slope, ran_3

      external
     $  ran_3
c
c H-mag distribution 
c
c Functions have been triple-checked as of 2018-05-04. OK.
      slope = h_params(3)
      h0s10 = 10.d0**(slope*h_params(1))
      h1s10 = 10.d0**(slope*h_params(2))
      random=ran_3(seed)
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
     $  xi1, ran_3

      external
     $  ran_3
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
      random=ran_3(seed)
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
     $  ran_3

      external
     $  ran_3
c
c H-mag distribution
c
c Functions have been triple-checked as of 2018-05-04. OK.
c
c In order to not have the ramp-up found in {\it size_dist_two}, I use
c an exponential with no lower limit on H. What sets the lower limit of
c the faint end is the selection of the correct slope at H_break (see
c test below). The test is based on the number of objects in each
c portion of the distribution.
c
c Since we want N(<H) = A 10^{\alpha_1 H} for H < H_b
c and           N(<H) = B 10^{\alpha_2 H} for H >= H_b
c and N(<H) continuous at H_b, we have N(<H_b) = A 10^{\alpha_1 H_b} =
c B 10^{\alpha_2 H_b}.
c
c The total number of objects is N(<H_max) = B 10^{\alpha_2 H_max}, and
c the number of objects bigger than H_b is N(<H_b) = B 10^{\alpha_2 H_b}. 
c So the fraction of objects bigger than H_b is 10^{\alpha_2 H_b} /
c 10^{\alpha_2 H_max} = h1s2 / h2s2 = h12s2.
c
c Therefore, if 'random' < h12s2, then we have a big object, and
c otherwise a small on. 'random' needs to be rescaled to go to 1.
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

      end

      real*8 function size_dist_n_straight (seed, h_params, nn)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This function draws a number according to an n-slope exponential
c cumulative distribution specified by "h_params" (n straight lines in
c semilog for cumulative):
c
c   P(<=h) = A_k 10^{(h*\alpha_k)}
c
c with H_{k-1} < h <= H_k,
c
c with
c
c \alpha_k = h_params(k+n)
c H_k = h_params(k)
c
c for k in [1; n ]. h_params(0) is not present and implicit at -\infty.
c This corresponds to dropping the lower limit as we want straight lines.
c
c The function is continuous at H_k = h_params(k) ofr all k's:
c   A_k 10^{H_k*\alpha_k} = A_{k+1} 10^{H_k*\alpha_{k+1}}
c
c for k in [1; n-1].
c
c To avoid using allocatable arrays and dynamical allocation, I restrict
c the number of slopes to be <= 10.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : September 2019
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     seed  : seed for the random number generator (I4)
c     h_params: parameters for the distribution (3*R8)
c     nn    : number of different slopes (I4)
c
c OUTPUT
c     size_dist_n_straight: Random value of H (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c Set of F2PY directives to create a Python module
c
Cf2py intent(in,out) seed
Cf2py intent(in) h_param
Cf2py intent(in) nn
c
      implicit none

      integer*4
     $  seed, n, nn, k

      real*8
     $  h_params(*), random, x(10), f(10), ran_3

      logical first

      save first

      external
     $  ran_3

      data first /.true./

      if (nn .gt. 10) then
         if (first) then
            first = .false.
            print *,
     $        'WARNING: number of breaks greater than allowed maximum'//
     $        ' 10, using only the 10 first values.'
         end if
         n = 10
      else
         n = nn
      end if
c
c H-mag distribution
c
c In order to not have the ramp-up found in {\it size_dist_two}, I use
c an exponential with no lower limit on H. What sets the lower limit of
c the faint end is the selection of the correct slope at H_k (see
c test below). The test is based on the number of objects in each
c portion of the distribution.
c
c The total number of objects is N(<H_max) = N(<H_n) = A_n 10^{H_n \alpha_n},
c and the number of objects bigger than the next break H_{n-1} is N(
c <H_{n-1}) = A_n 10^{H_{n-1} \alpha_n}. More generaly, the number of
c objects bigger than H_k is
c
c N_k = A_k 10^{H_k \alpha_k} = A_{k+1} 10^{H_k \alpha_{k+1}}
c
c Let's define the fraction of objects bigger than H_k in the objects
c bigger than H_{k+1} as:
c
c x_k = N_{k-1} / N_k = 10^{(H_{k-1}-H_k) \alpha_k}
c
c for k in [2; n ], and x_1 = 0.
c
c The fraction of object bigger than B_k compared to the total
c number of objects is
c
c f_k = N_k / N_n = (N_k / N_{k+1}) (N_{k+1} / N_{k+2}) ... (N_{n-1} / N_n)
c                 = x_{k+1} x_{k+2} ... x_n
c
c We define x_k as described above, then we set
c
c f_n = 1.
c f_{k-1} = x_k f_k, for k in [2; n ]
c
c Therefore, if (f_{k-1} < random <= f_k), then we have an object in
c range ]H_{k-1}; H_k]. 'random' needs to be rescaled to go to 1.
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
            size_dist_n_straight =
     $        log10(random/f(k)*10.d0**(h_params(k)*h_params(k+n)))
     $        /h_params(k+n)
            return
         end if
      end do

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

      real*8 FUNCTION ran_3(idum)
Cf2py intent(in,out) idum
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL*8 FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=abs(MSEED-abs(idum))
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
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
      END
