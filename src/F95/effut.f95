module effut

  use datadec

contains
  function eta_raw (eff_p, nr, mdum, rdum, maglim)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the efficiency at a given magnitude.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2004
! Version 2 : July 2004
!             Added 3 different types of efficiency functions, namely a
!             single and a double hyperbolic tangent, and a piecewise
!             linear function.
! Version 3 : June 2013
!             Added efficiency dependance on rate
! Version 4 : April 2014
!             Added limiting magnitude determination function of rate
! Version 5 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file)
! Ver f95-1 : March 2017
!             Identical to Version 5, but for F95
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     eff_p : Set of efficiency functions for rate ranges (n*eff_r)
!     nr    : Number of efficiency functions (I4)
!     mdum  : magnitude (R8)
!     rdum  : rate of motion (R8)
!
! OUTPUT
!     maglim: Limiting magnitude at given rate (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) eff_p
!f2py intent(in) nr
!f2py intent(in) mdum
!f2py intent(in) rdum
!f2py intent(out) maglim

    implicit none
    integer :: ilo, ihi, i, nr, ir
    type(t_eff_r) :: eff_p(*)
    real (kind=8) :: m, x, mytanh, mdum, rdum, r, maglim, eta_raw

    mytanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))

! Retrieve magnitude and rate of motion
    m = dmax1(0.d0, mdum)
    r = rdum

! Determine which efficiency function to use
    ir = 0
500 continue
      ir = ir + 1
      if ((r-eff_p(ir)%min)*(r-eff_p(ir)%max) .le. 0.d0) goto 600
      if (ir .lt. nr) goto 500
510 continue
    eta_raw = 0.d0
    maglim = 0.d0
    return
600 continue
! Report limiting magnitude for this rate range
    maglim = eff_p(ir)%mag_lim
!    write (18, *) 'Rate number: ', ir

! This is the direct piecewise linear function from a lokup table
    if (eff_p(ir)%n .gt. 0) then

! If off bins, then flat continuation.
       if (m .lt. eff_p(ir)%b(1)) then
          eta_raw = eff_p(ir)%e(1)
          return
       else if (m .gt. eff_p(ir)%b(eff_p(ir)%n)) then
          eta_raw = 0.d0
          return
       else

! Linear interpolation of table.
          ilo = 1
          ihi = eff_p(ir)%n
1000      continue
          if (ihi - ilo .gt. 1) then
             i = (ihi + ilo)/2
             if (eff_p(ir)%b(i) .lt. m) then
                ilo = i
             else if (eff_p(ir)%b(i) .gt. m) then
                ihi = i
             else
                eta_raw = eff_p(ir)%e(i)
                return
             end if
             goto 1000
          end if
          eta_raw = eff_p(ir)%e(ilo) + (eff_p(ir)%e(ihi) - eff_p(ir)%e(ilo))* &
               (m - eff_p(ir)%b(ilo))/(eff_p(ir)%b(ihi) - eff_p(ir)%b(ilo))
       end if

! This is a single hyperbolic tangent function.
! \begin{equation}
! (A/2) * (1. - tanh((R-R_c)/d))
! \end{equation}
    else if (eff_p(ir)%n .eq. -1) then
       eta_raw = eff_p(ir)%e(1)/2.d0 * &
            (1.d0 - mytanh((m - eff_p(ir)%e(2))/eff_p(ir)%e(3)))

! This is a double hyperbolic tangent function.
! \begin{equation}
! (A/4) * (1. - tanh((R-R_c)/d1)) * (1. - tanh((R-R_c)/d2))
! \end{equation}
    else if (eff_p(ir)%n .eq. -2) then
       eta_raw = eff_p(ir)%e(1)/4.d0 * &
            (1.d0 - mytanh((m - eff_p(ir)%e(2))/eff_p(ir)%e(3))) * &
            (1.d0 - mytanh((m - eff_p(ir)%e(2))/eff_p(ir)%e(4)))
!       write (18, *) eff_p(ir)%e(1), eff_p(ir)%e(2), eff_p(ir)%e(3), &
!            eff_p(ir)%e(4)
!       write (18, *) m
!       write (18, *) (m - eff_p(ir)%e(2))/eff_p(ir)%e(4), &
!            mytanh((m - eff_p(ir)%e(2))/eff_p(ir)%e(4)), &
!            exp((m - eff_p(ir)%e(2))/eff_p(ir)%e(4))

! This is a piecewize linear function.
! \begin{eqnarray}
! A & {\rm if} & m < R_1 \\
! \frac{(m - R_2) A}{R_1 - R_2} & {\rm if} & R_1 \le m < R_2 \\
! 0 & {\rm if} & m \ge R_2
! \end{eqnarray}
    else if (eff_p(ir)%n .eq. -3) then
       if (m .lt. eff_p(ir)%e(2)) then
          eta_raw = eff_p(ir)%e(1)
       else if (m .lt. eff_p(ir)%e(3)) then
          eta_raw = (m - eff_p(ir)%e(3))*eff_p(ir)%e(1) / &
               (eff_p(ir)%e(2) - eff_p(ir)%e(3))
       else
          eta_raw = 0.d0
       end if

! This is the SKADS defined function
! \begin{equation}
! (A - c * (R-21.0)^2) / (1. + exph((R-R_c)/d))
! \end{equation}
    else if (eff_p(ir)%n .eq. -4) then
       if (m .lt. 21.d0) then
          eta_raw = eff_p(ir)%e(1)
       else
          eta_raw = (eff_p(ir)%e(1) - eff_p(ir)%e(2) * (m - 21.d0)**2) / &
               (1.d0 + exp((m - eff_p(ir)%e(3))/eff_p(ir)%e(4)))
       end if

! Unsupported efficiency function type.
    else
       write (6, *) 'Got efficiency function type ', eff_p(ir)%n
       write (6, *) 'Should be >0, -1, -2 or -3.'
       stop 'Something is wrong with this. Aborting.'
    end if

    return
  end function eta_raw

  function eta (eff_p, nr, mdum, rdum, maglim)

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the efficiency at a given magnitude, keeping
! only the part that can be trusted for actual detectability of the
! theoretical magnitude ($\eta > 0.4$).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
! Version 3 : June 2013
!             Added efficiency dependance on rate
! Version 4 : April 2014
!             Added limiting magnitude determination function of rate
! Version 5 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     eff_p : Set of efficiency functions for rate ranges (n*eff_r)
!     nr    : Number of efficiency fucntions (I4)
!     mdum  : magnitude (R8)
!     rdum  : rate of motion (R8)
!
! OUTPUT
!     maglim: Limiting magnitude at given rate (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) eff_p
!f2py intent(in) nr
!f2py intent(in) mdum
!f2py intent(in) rdum
!f2py intent(out) maglim

    implicit none
    integer :: eff_n(n_r_max), nr
    type(t_eff_r) :: eff_p(*)
    real (kind=8) :: mdum, rdum, maglim, eta
    real (kind=8), save :: lim

    data lim /0.01d0/

    eta = eta_raw(eff_p, nr, mdum, rdum, maglim)
    if (eta .lt. lim) eta = 0.d0

    return
  end function eta

end module effut
