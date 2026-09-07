program test_eta_rate
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Unit tests for rate-dependent m0 in the single hyperbolic-tangent
! efficiency function:
!   eta = (eta0/2) * (1 - tanh((m - m0)/w))
!   m0  = m00 + alpha_rate * rate_asphr
! where rate_asphr is the total on-sky rate in arcsec/hour.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  use parameters
  use datadec
  use effut
  use getsur

  implicit none

  integer, parameter :: lun_in = 21
  type(t_eff_r), dimension(1) :: eff_p
  type(t_charact) :: c
  integer :: ierr, nfail
  real (kind=8) :: eta0, m00, w, alpha, mag, maglim
  real (kind=8) :: rate1_asp, rate2_asp, rate1, rate2
  real (kind=8) :: eta1, eta2, exp1, exp2, eta_legacy
  real (kind=8) :: m0_1, m0_2

  nfail = 0

  eta0 = 0.90d0
  m00 = 24.00d0
  w = 0.50d0
  alpha = 0.10d0
  mag = 24.00d0
  rate1_asp = 1.0d0
  rate2_asp = 3.0d0
  rate1 = rate1_asp*24.d0/3600.d0*drad
  rate2 = rate2_asp*24.d0/3600.d0*drad

  ! ---- In-memory evaluation with alpha_rate ----
  eff_p(1)%min = 0.d0
  eff_p(1)%max = 20.d0*24.d0/3600.d0*drad
  eff_p(1)%mag_lim = 24.5d0
  eff_p(1)%n = -1
  eff_p(1)%e = 0.d0
  eff_p(1)%e(1) = eta0
  eff_p(1)%e(2) = m00
  eff_p(1)%e(3) = w
  eff_p(1)%e(4) = alpha

  eta1 = eta_raw(eff_p, 1, mag, rate1, maglim)
  eta2 = eta_raw(eff_p, 1, mag, rate2, maglim)

  m0_1 = m00 + alpha*rate1_asp
  m0_2 = m00 + alpha*rate2_asp
  exp1 = eta0/2.d0*(1.d0 - dtanh((mag - m0_1)/w))
  exp2 = eta0/2.d0*(1.d0 - dtanh((mag - m0_2)/w))

  call check('in-memory eta at 1"/hr', eta1, exp1, 1.d-12, nfail)
  call check('in-memory eta at 3"/hr', eta2, exp2, 1.d-12, nfail)
  if (eta2 .le. eta1) then
     write (6, *) 'FAIL: faster rate should raise m0 and increase eta at fixed mag'
     nfail = nfail + 1
  else
     write (6, *) 'PASS: faster rate increases eta at fixed mag (brighter m0)'
  end if

  ! ---- alpha_rate = 0 recovers constant-m0 behaviour ----
  eff_p(1)%e(4) = 0.d0
  eta_legacy = eta_raw(eff_p, 1, mag, rate2, maglim)
  exp1 = eta0/2.d0*(1.d0 - dtanh((mag - m00)/w))
  call check('alpha=0 ignores rate', eta_legacy, exp1, 1.d-12, nfail)
  eta1 = eta_raw(eff_p, 1, mag, rate1, maglim)
  call check('alpha=0 same at both rates', eta1, eta_legacy, 1.d-12, nfail)

  ! ---- read_eff: 4-parameter single_param ----
  call read_eff('data/rate_m0.eff', lun_in, c, ierr)
  if (ierr .ne. 0) then
     write (6, *) 'FAIL: read_eff rate_m0.eff ierr=', ierr
     nfail = nfail + 1
  else
     call check('read alpha_rate', c%eff_p(1)%e(4), alpha, 1.d-12, nfail)
     call check('read eta0', c%eff_p(1)%e(1), eta0, 1.d-12, nfail)
     call check('read m00', c%eff_p(1)%e(2), m00, 1.d-12, nfail)
     call check('read w', c%eff_p(1)%e(3), w, 1.d-12, nfail)
     eta1 = eta_raw(c%eff_p, c%nr, mag, rate1, maglim)
     eta2 = eta_raw(c%eff_p, c%nr, mag, rate2, maglim)
     call check('file-backed eta at 1"/hr', eta1, &
          eta0/2.d0*(1.d0 - dtanh((mag - (m00+alpha*rate1_asp))/w)), &
          1.d-12, nfail)
     call check('file-backed eta at 3"/hr', eta2, &
          eta0/2.d0*(1.d0 - dtanh((mag - (m00+alpha*rate2_asp))/w)), &
          1.d-12, nfail)
  end if

  ! ---- read_eff: legacy 3-parameter single_param ----
  call read_eff('data/rate_m0_legacy.eff', lun_in, c, ierr)
  if (ierr .ne. 0) then
     write (6, *) 'FAIL: read_eff rate_m0_legacy.eff ierr=', ierr
     nfail = nfail + 1
  else
     call check('legacy alpha defaults to 0', c%eff_p(1)%e(4), 0.d0, &
          1.d-15, nfail)
     eta1 = eta_raw(c%eff_p, c%nr, mag, rate1, maglim)
     eta2 = eta_raw(c%eff_p, c%nr, mag, rate2, maglim)
     call check('legacy eta independent of rate (1)', eta1, exp1, 1.d-12, nfail)
     call check('legacy eta independent of rate (3)', eta2, exp1, 1.d-12, nfail)
  end if

  if (nfail .eq. 0) then
     write (6, *) 'ALL TESTS PASSED'
     stop 0
  else
     write (6, *) 'FAILURES: ', nfail
     stop 1
  end if

contains

  subroutine check(label, got, expected, tol, nfail)
    character(*), intent(in) :: label
    real (kind=8), intent(in) :: got, expected, tol
    integer, intent(inout) :: nfail
    if (abs(got - expected) .gt. tol) then
       write (6, *) 'FAIL:', label, 'got', got, 'expected', expected
       nfail = nfail + 1
    else
       write (6, *) 'PASS:', label
    end if
  end subroutine check

end program test_eta_rate
