"""
Functions us in analysis of OSSOS H distribution.

"""
import numpy
from scipy import stats
ln10 = numpy.log(10.)


def broken_plaw(H, a1, a2, Ho, Hb):
    """
    LF from Fraser et al. 2014
    :param H: h magnitude
    :param a1: log slope bright of break)
    :param a2: log slope faint of break)
    :param Ho: normalization of bright component
    :param Hb: normalization of faint component
    :return: N

    As presented here this the nominal cummulative form of the function.
    """
    N = 10**(a1*(H-Ho))
    N[H > Hb] = 10**(a2*(H[H > Hb]-Ho)+(a1-a2)*(Hb-Ho))
    return N


def HBO(alpha_SI, beta_SI, A, B):
    """
    Convert the 'A' and 'B' parameters in variably tappered funciton to H1 and H2.
    :param alpha_SI:
    :param beta_SI:
    :param A:
    :param B:
    :return: H1, H2

    The functional form used in the analysis users constants 'A' and 'B' rather than 10**H1 and 10**H2 as the fits are more well
    behaved.  This function converts from A/B to the H normalization.
    """
    return -5/(3*alpha_SI) * numpy.log10(A), 5/(3*beta_SI) * numpy.log10(B)


def variably_tapered(h, A, B, alpha_SI, beta_SI):
    """
    Exponentially tappered powerlaw size distribution re-expressed in H magnitude space.
    :param A: The normalization of the asymtoptic exponential
    :param B: The point where the exponential taper begins.
    :param alpah_SI: exponential asymptotic value.
    :param beta_SI: exponent of the taper.

    Derived from Schmit et al.  This form of the LF mimics behaviour seen in the streaming instability planetesimal formation process.
    """
    return A*10.**(alpha_SI*3.*h/5.)*numpy.exp(-B*10.**(-beta_SI*3.*h/5.))


def variably_tapered2(h, Ho, Hb, alpha_SI, beta_SI):
    """
    Exponentially tapered exponential size distribution re-expressed in H magnitude space.
    """
    return 10.**(alpha_SI*3./5.*(h-Ho))*numpy.exp(-10.**(-beta_SI*3./5*(h-Hb)))


def variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI):
    """
    Differential form of the variably_tapere
    :param h:
    :param Dp:
    :param E:
    :param alpha_SI:
    :param beta_SI:
    :return:
    """

    return (3./5.*ln10)*Dp*10.**(alpha_SI*3.*h/5.)*(alpha_SI+beta_SI*E*10.**(-beta_SI*3.*h/5.))*numpy.exp(-E*10.**(-beta_SI*3.*h/5.))


def log_variably_tapered(h, Dp, E, alpha_SI, beta_SI):
    return numpy.log10(variably_tapered(h, Dp, E, alpha_SI, beta_SI))


def log_variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI):
    return numpy.log10(variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI))


def likelihood(H, bias, A=0.88, B=450.0,  alpha_SI=0.67, beta_SI=0.65):
    """

    :param H: Vector of measured H values
    :param bias: Vector of bias on detection of a given object


    :return:
    """
    x = numpy.linspace(H.min()-3.0, H.max())
    dx = x[1]-x[0]
    cdf = variably_tapered(x, A, B, alpha_SI, beta_SI)
    pdf = numpy.diff(cdf)/dx
    ht = numpy.linspace(H.min(), H.max())
    dh = ht[1]-ht[0]
    dy = bias*numpy.interp(H, x[1:]-dx/2.0, pdf)
    dN = numpy.interp(ht, H, dy)
    if False:
        plt.clf()
        plt.plot(H, numpy.interp(H, x, cdf), '-')
        plt.plot(x[1:], (dx*pdf).cumsum(), ':k')
        plt.plot(ht, dh*dN.cumsum(), ':k')
        plt.plot(H, numpy.ones(len(H)).cumsum(), '-')
        plt.yscale('log')
        plt.ylim(1, 10000)
        plt.show()

    N = dh*dN.sum()
    l = -N + numpy.sum(numpy.log(dy/len(H)))
    return dy, N, l


def likelihood2(H, bias, Ho=0, Hb=7,  alpha_SI=0.67, beta_SI=0.65):
    """

    :param H: Vector of measured H values
    :param bias: Vector of bias on detection of a given object


    :return:
    """
    x = numpy.linspace(H.min()-3.0, H.max()+3.0)
    dx = x[1]-x[0]
    cdf = variably_tapered2(x, Ho, Hb, alpha_SI, beta_SI)
    pdf = numpy.diff(cdf)/dx
    ht = numpy.linspace(H.min(), H.max())
    dh = ht[1]-ht[0]
    dy = bias*numpy.interp(H, x[1:]-dx/2.0, pdf)
    dN = numpy.interp(ht, H, dy)

    N = dh*dN.sum()
    l = -N + numpy.sum(numpy.log(dy/len(H)))
    return dy, N, l

def fraser(H, a1, a2, Ho, Hb, dh):
    """
    LF from Fraser et al. 2014
    :param H: h magnitude
    :param a1: log slope bright of break)
    :param a2: log slope faint of break)
    :param Ho: normalization of bright component
    :param Hb: normalization of faint component
    :return: N
    """
    # print(f"alpha1: {a1}, alpha2: {a2}, Ho: {Ho}, Hb:{Hb}")
    # N[H > Hb] = a2*ln10*10**(a2*(H[H > Hb]-Ho)+(a1-a2)*(Hb-Ho))
    # 10**(a2*H-a2*Ho+a1*Hb-a1*Ho-a2*Hb+a2*Ho)
    # 10**(a2*(H-Hb)+a1*(Hb-Ho))
    # 10**(a1*(Hb-Ho)) * 10**(a2*(H-Hb))
    N = dh*a1*ln10*10**(a1*(H-Ho))
    C = 10**(a1*(Hb-Ho))
    N[H >= Hb] = dh*a1*ln10*C*10**(a2*(H[H > Hb]-Hb))
    print(N[H<Hb][-3:], N[H>=Hb][0])
    return N


def double_plaw(H, simga_23=0.68, a1=1.36, a2=0.38, R_eq=22.8, dh=0.1):
    """
    Double powerlaw form from B14
    :param H:
    :param simga_23:
    :param a1:
    :param a2:
    :return: simga(H)
    """
    r = H + 10*numpy.log10(42.0)
    C = 10**((a2-a2)*(R_eq - 23))
    return dh*((1+C)*simga_23/(10**(-a1*(r-23)) + C*10**(-a2*(r-23))))


def rolling_plaw(H, sigma_23, a1, a2):
    """
    Rolling power law form form Bernstein et al. 2004
    :param H: numpy array of H mags.
    :param sigma_23: normalization at R=23
    :param a1: bright slope.
    :param a2: faint slope.
    :return: N
    """
    r = H + 10*numpy.log10(44)
    return sigma_23 * 10**(a1*(r-23)+a2*(r-23)**2)




def ll(params, H, bias, limits):
    """
    Compute the log-likelihood [setup to work with eemcc
    :param params: function parameters
    :param H: numpy array of H magnitudes
    :param bias: numpy array of detection bias
    :param limits: bounds on params, ll = -np.inf outside this range.
    :return: log of likelihood
    """
    for idx in ['Ho', 'Hb', 'alpha_SI', 'beta_SI']:
        if not limits[idx][0] < params[idx] < limits[idx][1]:
            return -numpy.inf
    return likelihood2(H, bias, alpha_SI=params[2], Hb=params[1], Ho=params[0], beta_SI=params[3])[2]


def poisson_range(k, prob):
    """
    Given an measured count rate, k, and an expectation of a poisson process return estimates of 'mu' that
    are consistent with 'k' being measured within probability prob.

    i.e. if we measure k objects compute a bunch of mu and return the range of mu that are consistent within prob.
    """
    # mu holds the  range of plausible estimates for actual mu
    mu = numpy.arange(max(0, k.min()-20*(k.min())**0.5), k.max()+20*(k.max())**0.5+10, 10+20*k.max()**0.5/100.)
    upper = stats.poisson.ppf(0.5-prob/2., mu)
    lower = stats.poisson.ppf(0.5+prob/2., mu)
    return numpy.interp(k, lower, mu), numpy.interp(k, upper, mu)


if __name__ == '__main__':
    print("This is the result from the OSSOS SFD papers.")
    from matplotlib import pyplot as plt
    Ho = -2.6
    Hb = 8.1 
    beta_SI = 0.42
    alpha_SI = 0.677
    x = numpy.arange(4.5, 17, 0.1)
    # x = numpy.array([8.66,12,17.0])
    # From Kavelaars et al. (2021)
    # alpha = 3/5 alpha_SI => alpha_SI = 5*alpha/3 
    alpha_SI = 5*0.4/3
    beta_SI = 0.42
    Ho = -2.6
    Hb = 8.1
    cdf = variably_tapered2(x, Ho, Hb, alpha_SI, beta_SI)
    plt.plot(x, cdf)
    plt.xlabel('H (r)')
    plt.ylabel('N(H<r)')
    plt.yscale('log')
    plt.savefig('SFD.pdf')
