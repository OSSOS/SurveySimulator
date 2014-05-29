      real*8 function eta_raw (nb_max, rates, nr, eff_n, eff_b,
     $  eff_m, mdum, rdum, ml, maglim)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the efficiency at a given magnitude.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : April 2004
c Version 2 : July 2004
c             Added 3 different types of efficiency functions, namely a
c             single and a double hyperbolic tangent, and a piecewise
c             linear function.
c Version 3 : June 2013
c             Added efficiency dependance on rate
c Version 4 : April 2014
c             Added limiting magnitude determination function of rate
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nb_max: Maximum number of bin in efficiency function (I4)
c     rates : Rates limits for efficiency function (2,n*R8)
c     nr    : Number of efficiency fucntions (I4)
c     eff_n : Number of efficiency bins (n*I4)
c     eff_b : Magnitude bin center (nb_max,n*R8)
c     eff_m : Efficiency for that magnitude (nb_max,n*R8)
c     mdum  : magnitude (R8)
c     rdum  : rate of motion (R8)
c     ml    : Limiting magnitudes for rate ranges (n*R8)
c
c OUTPUT
c     maglim: Limiting magnitude at given rate (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      integer*4
     $  ilo, ihi, i, eff_n(*), nr, nb_max, ir

      real*8
     $  m, x, tanh, eff_b(nb_max,*), eff_m(nb_max,*), mdum, rdum,
     $  rates(2,*), r, ml(*), maglim

      tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))

c Retrieve magnitude and rate of motion
      m = dmax1(0.d0, mdum)
      r = rdum

c Determine which efficiency function to use
      ir = 0
 500  continue
         ir = ir + 1
         if ((r-rates(1,ir))*(r-rates(2,ir)) .le. 0.d0) goto 600
         if (ir .lt. nr) goto 500
 510  continue
      eta_raw = 0.d0
      maglim = 0.d0
      return
 600  continue
c Report limiting magnitude for this rate range
      maglim = ml(ir)

c This is the direct piecewise linear function from a lokup table
      if (eff_n(ir) .gt. 0) then

c If off bins, then flat continuation.
         if (m .lt. eff_b(1,ir)) then
            eta_raw = eff_m(1,ir)
            return
         else if (m .gt. eff_b(eff_n(ir),ir)) then
            eta_raw = 0.d0
            return
         else

c Linear interpolation of table.
            ilo = 1
            ihi = eff_n(ir)
 1000       continue
            if (ihi - ilo .gt. 1) then
               i = (ihi + ilo)/2
               if (eff_b(i,ir) .lt. m) then
                  ilo = i
               else if (eff_b(i,ir) .gt. m) then
                  ihi = i
               else
                  eta_raw = eff_m(i,ir)
                  return
               end if
               goto 1000
            end if
            eta_raw = eff_m(ilo,ir) + (eff_m(ihi,ir) - eff_m(ilo,ir))*
     $        (m - eff_b(ilo,ir))/(eff_b(ihi,ir) - eff_b(ilo,ir))
         end if

c This is a single hyperbolic tangent function.
c \begin{equation}
c (A/2) * (1. - tanh((R-R_c)/d))
c \end{equation}
      else if (eff_n(ir) .eq. -1) then
         eta_raw = eff_m(1,ir)/2.d0 *
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(3,ir)))

c This is a double hyperbolic tangent function.
c \begin{equation}
c (A/4) * (1. - tanh((R-R_c)/d1)) * (1. - tanh((R-R_c)/d2))
c \end{equation}
      else if (eff_n(ir) .eq. -2) then
         eta_raw = eff_m(1,ir)/4.d0 *
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(3,ir))) * 
     $     (1.d0 - tanh((m - eff_m(2,ir))/eff_m(4,ir)))

c This is a piecewize linear function.
c \begin{eqnarray}
c A & {\rm if} & m < R_1 \\
c \frac{(m - R_2) A}{R_1 - R_2} & {\rm if} & R_1 \le m < R_2 \\
c 0 & {\rm if} & m \ge R_2
c \end{eqnarray}
      else if (eff_n(ir) .eq. -3) then
         if (m .lt. eff_m(2,ir)) then
            eta_raw = eff_m(1,ir)
         else if (m .lt. eff_m(3,ir)) then
            eta_raw = (m - eff_m(3,ir))*eff_m(1,ir)/
     $        (eff_m(2,ir) - eff_m(3,ir))
         else
            eta_raw = 0.d0
         end if

c This is the SKADS defined function
c \begin{equation}
c (A - c * (R-21.0)^2) / (1. + exph((R-R_c)/d))
c \end{equation}
      else if (eff_n(ir) .eq. -4) then
         if (m .lt. 21.d0) then
            eta_raw = eff_m(1,ir)
         else
            eta_raw = (eff_m(1,ir) -
     $        eff_m(2,ir) * (m - 21.d0)**2) /
     $        (1.d0 + exp((m - eff_m(3,ir))/eff_m(4,ir)))
         end if

c Unsupported efficiency function type.
      else
         write (6, *) 'Got efficiency function type ', eff_n(ir)
         write (6, *) 'Should be >0, -1, -2 or -3.'
         stop 'Something is wrong with this. Aborting.'
      end if

      return
      end

      real*8 function eta_trust (nb_max, rates, nr, eff_n, eff_b,
     $  eff_m, mdum, rdum, ml, maglim)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the efficiency at a given magnitude, keeping
c only the part that can be trusted for actual detectability of the
c theoretical magnitude ($\eta > 0.2$).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c Version 3 : June 2013
c             Added efficiency dependance on rate
c Version 4 : April 2014
c             Added limiting magnitude determination function of rate
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nb_max: Maximum number of bin in efficiency function (I4)
c     rates : Rates limits for efficiency function (2,n*R8)
c     nr    : Number of efficiency fucntions (I4)
c     eff_n : Number of efficiency bins (n*I4)
c     eff_b : Magnitude bin center (nb_max,n*R8)
c     eff_m : Efficiency for that magnitude (nb_max,n*R8)
c     mdum  : magnitude (R8)
c     rdum  : rate of motion (R8)
c     ml    : Limiting magnitudes for rate ranges (n*R8)
c
c OUTPUT
c     maglim: Limiting magnitude at given rate (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      integer*4
     $  eff_n(*), nr, nb_max

      real*8
     $  eta_raw, eff_b(nb_max,*), eff_m(nb_max,*), mdum, rdum, limtrust,
     $  rates(2,*), ml(*), maglim

      external
     $  eta_raw

      data
     $  limtrust /0.2d0/

      save limtrust

      eta_trust = eta_raw(nb_max, rates, nr, eff_n, eff_b, eff_m,
     $  mdum, rdum, ml, maglim)
      if (eta_trust .lt. limtrust) eta_trust = 0.d0

      return
      end

      real*8 function eta (nb_max, rates, nr, eff_n, eff_b,
     $  eff_m, mdum, rdum, ml, maglim)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine computes the efficiency at a given magnitude, keeping
c only the part that can be trusted for actual detectability of the
c theoretical magnitude ($\eta > 0.4$).
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : October 2006
c Version 3 : June 2013
c             Added efficiency dependance on rate
c Version 4 : April 2014
c             Added limiting magnitude determination function of rate
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nb_max: Maximum number of bin in efficiency function (I4)
c     rates : Rates limits for efficiency function (2,n*R8)
c     nr    : Number of efficiency fucntions (I4)
c     eff_n : Number of efficiency bins (n*I4)
c     eff_b : Magnitude bin center (nb_max,n*R8)
c     eff_m : Efficiency for that magnitude (nb_max,n*R8)
c     mdum  : magnitude (R8)
c     rdum  : rate of motion (R8)
c     ml    : Limiting magnitudes for rate ranges (n*R8)
c
c OUTPUT
c     maglim: Limiting magnitude at given rate (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

      implicit none

      integer*4
     $  eff_n(*), nr, nb_max

      real*8
     $  eta_raw, eff_b(nb_max,*), eff_m(nb_max,*), mdum, rdum, lim,
     $  rates(2,*), ml(*), maglim

      external
     $  eta_raw

      data
     $  lim /0.01d0/

      save lim

      eta = eta_raw(nb_max, rates, nr, eff_n, eff_b, eff_m, mdum, rdum,
     $  ml, maglim)
      if (eta .lt. lim) eta = 0.d0

      return
      end
