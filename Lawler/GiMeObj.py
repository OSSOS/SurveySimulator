import numpy as np

deg2rad = np.pi / 180.0  # Converts degrees into radians.

rnd_gen = None  # Creates empty variable to hold the random number generator once the seed in set.

# Creates empty global variables to house the object class instances generated in initialize().
obj1 = None
obj2 = None
obj3 = None
obj4 = None
obj5 = None

alpha_f, alpha_b, hmax, hmin, hbreak, contrast = \
            np.genfromtxt('H_dist_vars.txt', unpack=True)


class Obj:
    """
    Class used to create and store the objects generated and passed by the GiMeObj module into the main Driver.py
    module that executes the survey simulator code.
    """

    def __init__(self, size=10 ** 6):
        """
        Initializes an object to contain the parameter arrays based on called distributions.

        Any variables that will be the same for all object populations can be set here instead of being assigned in
        the individual draw methods.

        Parameters
        ----------
        size = (optional) Determines size of the arrays to be generated (default=10^6).
        """

        self.size = size

        self.index = 0  # Initializes an index used to count over the elements of the arrays.

        self.comments = None

        self.a = None
        self.e = None
        self.i = None
        self.node = None
        self.peri = None
        self.m = None
        self.epoch = None
        self.h = None
        self.color = None
        self.gb = None
        self.ph = None
        self.period = None
        self.amp = None
        self.resamp = None

    def draw_symm(self):
        """
        Initializes an Obj class object used to contain the arrays of a given size filled with the population
        distributions of the various orbital parameters for a single orbital population (e.g. the plutinos).

        This distribution represents a toy model of objects in the symmetric resonance island with Neptune.
        """

        self.comments = "Symmetric"

        self.a = self.uniform(47.8, 48.8)
        # self.e = self.gauss_range(0.18, 0.06, 0, 1)
        # self.i = self.sin_gauss_range(0 * deg2rad, 16 * deg2rad, 0 * deg2rad, 90 * deg2rad)
        self.e = self.constant(1-30.1/47.8)
        self.i = self.constant(0.001)
        self.node = self.uniform(0, 2*np.pi)  # random value between 0 and 2pi
        self.m = self.uniform(0, 2*np.pi)   # random value between 0 and 2pi
        self.epoch = self.constant(2456293.5)
        self.h = self.power_knee_divot(alpha_b, hmax, hmin, hbreak, alpha_f, contrast)
        self.gb = self.constant(-0.12)  # Opposition surge effect
        self.ph = self.constant(0.00)  # Initial phase of lightcurve
        self.period = self.constant(0.60)  # Period of lightcurve
        self.amp = self.constant(0.00)  # Amplitude of lightcurve (peak-to-peak)
        self.resamp = self.linear_peak(0, 5, 10)

        # The following section of code builds the 2-dimensional color array of dimension (size, 10)
        #
        # This color array is set up for the r band, as such r-x is necessarily 0.
        #
        # color : Array of colors (10*R8)
        #     objX.color[index,0] : g-x
        #     objX.color[index,1] : r-x
        #     objX.color[index,2] : i-x
        #     objX.color[index,3] : z-x
        #     objX.color[index,4] : u-x
        #     objX.color[index,5] : V-x
        #     objX.color[index,6] : B-x
        #     objX.color[index,7] : R-x
        #     objX.color[index,8] : I-x
        #     objX.color[index,9] : extra
        #         (This tenth "extra" parameter currently has no physical meaning, but in included because the fortran
        #         code expects an array of this size. this parameter is intended for future expansions.)

        g_x = self.gauss(0.7, 0.2)  # g-r color sampled from gaussian with mu 0.7 and std 0.2 (from CFEPS data)
        r_x = self.constant(0)
        i_x = self.constant(0)
        z_x = self.constant(0)
        u_x = self.constant(0)
        V_x = self.constant(0)
        B_x = self.constant(0)
        R_x = self.constant(-0.1)
        I_x = self.constant(0)
        extra = self.constant(0)

        self.color = np.column_stack((g_x, r_x, i_x, z_x, u_x, V_x, B_x, R_x, I_x, extra))

        # choose value of phi sinusoidally from resamp
        phi = np.pi + np.sin(self.uniform(0, 2*np.pi))*self.resamp*deg2rad
        lambda_n = 333.5078 * deg2rad  # Neptune's mean longitude on 1 Jan 2013

        self.peri = (1 / 1 * (phi - 2 * self.m) - self.node + lambda_n) % (2 * np.pi)

    def draw_asymm_lead(self):
        """
        Initializes an Obj class object used to contain the arrays of a given size filled with the population
        distributions of the various orbital parameters for a single orbital population (e.g. the plutinos).

        This distribution represents a toy model of objects in the asymmetric leading resonance island with Neptune.
        """

        self.comments = "Asymm_lead"

        self.a = self.uniform(47.8, 48.8)
        # self.e = self.gauss_range(0.18, 0.06, 0, 1)
        # self.i = self.sin_gauss_range(0 * deg2rad, 16 * deg2rad, 0 * deg2rad, 90 * deg2rad)
        self.e = self.constant(1-30.1/47.8)
        self.i = self.constant(0.001)
        self.node = self.uniform(0, 2*np.pi)  # random value between 0 and 2pi
        self.m = self.uniform(0, 2*np.pi)   # random value between 0 and 2pi
        self.epoch = self.constant(2456293.5)
        self.h = self.power_knee_divot(alpha_b, hmax, hmin, hbreak, alpha_f, contrast)
        self.gb = self.constant(-0.12)  # Opposition surge effect
        self.ph = self.constant(0.00)  # Initial phase of lightcurve
        self.period = self.constant(0.60)  # Period of lightcurve
        self.amp = self.constant(0.00)  # Amplitude of lightcurve (peak-to-peak)
        self.resamp = self.linear_peak(0, 5, 10)

        # The following section of code builds the 2-dimensional color array of dimension (size, 10)
        #
        # This color array is set up for the r band, as such r-x is necessarily 0.
        #
        # color : Array of colors (10*R8)
        #     objX.color[index,0] : g-x
        #     objX.color[index,1] : r-x
        #     objX.color[index,2] : i-x
        #     objX.color[index,3] : z-x
        #     objX.color[index,4] : u-x
        #     objX.color[index,5] : V-x
        #     objX.color[index,6] : B-x
        #     objX.color[index,7] : R-x
        #     objX.color[index,8] : I-x
        #     objX.color[index,9] : extra
        #         (This tenth "extra" parameter currently has no physical meaning, but in included because the fortran
        #         code expects an array of this size. this parameter is intended for future expansions.)

        g_x = self.gauss(0.7, 0.2)  # g-r color sampled from gaussian with mu 0.7 and std 0.2 (from CFEPS data)
        r_x = self.constant(0)
        i_x = self.constant(0)
        z_x = self.constant(0)
        u_x = self.constant(0)
        V_x = self.constant(0)
        B_x = self.constant(0)
        R_x = self.constant(-0.1)
        I_x = self.constant(0)
        extra = self.constant(0)

        self.color = np.column_stack((g_x, r_x, i_x, z_x, u_x, V_x, B_x, R_x, I_x, extra))

        # choose value of phi sinusoidally from resamp
        phi = 80*deg2rad + np.sin(self.uniform(0, 2*np.pi))*self.resamp*deg2rad
        lambda_n = 333.5078 * deg2rad  # Neptune's mean longitude on 1 Jan 2013

        self.peri = (1 / 1 * (phi - 2 * self.m) - self.node + lambda_n) % (2 * np.pi)

    def draw_asymm_trail(self):
        """
        Initializes an Obj class object used to contain the arrays of a given size filled with the population
        distributions of the various orbital parameters for a single orbital population (e.g. the plutinos).

        This distribution represents a toy model of objects in the asymmetric trailing resonance island with Neptune.
        """

        self.comments = "Asymm_trail"

        self.a = self.uniform(47.8, 48.8)
        # self.e = self.gauss_range(0.18, 0.06, 0, 1)
        # self.i = self.sin_gauss_range(0 * deg2rad, 16 * deg2rad, 0 * deg2rad, 90 * deg2rad)
        self.e = self.constant(1-30.1/47.8)
        self.i = self.constant(0.001)
        self.node = self.uniform(0, 2*np.pi)  # random value between 0 and 2pi
        self.m = self.uniform(0, 2*np.pi)   # random value between 0 and 2pi
        self.epoch = self.constant(2456293.5)
        self.h = self.power_knee_divot(alpha_b, hmax, hmin, hbreak, alpha_f, contrast)
        self.gb = self.constant(-0.12)  # Opposition surge effect
        self.ph = self.constant(0.00)  # Initial phase of lightcurve
        self.period = self.constant(0.60)  # Period of lightcurve
        self.amp = self.constant(0.00)  # Amplitude of lightcurve (peak-to-peak)
        self.resamp = self.linear_peak(0, 5, 10)

        # The following section of code builds the 2-dimensional color array of dimension (size, 10)
        #
        # This color array is set up for the r band, as such r-x is necessarily 0.
        #
        # color : Array of colors (10*R8)
        #     objX.color[index,0] : g-x
        #     objX.color[index,1] : r-x
        #     objX.color[index,2] : i-x
        #     objX.color[index,3] : z-x
        #     objX.color[index,4] : u-x
        #     objX.color[index,5] : V-x
        #     objX.color[index,6] : B-x
        #     objX.color[index,7] : R-x
        #     objX.color[index,8] : I-x
        #     objX.color[index,9] : extra
        #         (This tenth "extra" parameter currently has no physical meaning, but in included because the fortran
        #         code expects an array of this size. this parameter is intended for future expansions.)

        g_x = self.gauss(0.7, 0.2)  # g-r color sampled from gaussian with mu 0.7 and std 0.2 (from CFEPS data)
        r_x = self.constant(0)
        i_x = self.constant(0)
        z_x = self.constant(0)
        u_x = self.constant(0)
        V_x = self.constant(0)
        B_x = self.constant(0)
        R_x = self.constant(-0.1)
        I_x = self.constant(0)
        extra = self.constant(0)

        self.color = np.column_stack((g_x, r_x, i_x, z_x, u_x, V_x, B_x, R_x, I_x, extra))

        # choose value of phi sinusoidally from resamp
        phi = 280*deg2rad + np.sin(self.uniform(0, 2*np.pi))*self.resamp*deg2rad
        lambda_n = 333.5078 * deg2rad  # Neptune's mean longitude on 1 Jan 2013

        self.peri = (1 / 1 * (phi - 2 * self.m) - self.node + lambda_n) % (2 * np.pi)

    def constant(self, constant):
        """
        Obj method that creates an array of the appropriate size filled with a constant value.

        Parameters
        ----------
        constant = The chosen constant value.
        """

        return_array = np.full(self.size, constant)

        return return_array

    def uniform(self, minimum=0., maximum=1.):
        """
        Obj method that creates an array of the appropriate size filled with random variables on the interval [0,1).

        Parameters
        ----------
        minimum = (optional) Minimum value to be randomly generated (default=0).
        maximum  = (optional) Maximum value to be randomly generated (default=1).
        """

        return_array = rnd_gen.uniform(minimum, maximum, self.size)

        return return_array

    def gauss(self, mu, sigma):
        """
        Obj method that creates an array of the appropriate size filled with random variables in a gaussian
        distribution.

        Parameters
        ----------
        mu = The mean of the gaussian distribution.
        sigma = The standard deviation of the gaussian distribution.
        """

        return_array = rnd_gen.normal(mu, sigma, self.size)

        return return_array

    def gauss_range(self, mu, sigma, minimum, maximum):
        """
        Obj method used to create an array corresponding to a probability density of a gaussian distribution.
        (i.e. y=exp(-(x-mu)^2/(2*sigma^2))

        Parameters
        ----------
        mu = The mean of the gaussian distribution.
        sigma = The standard deviation of the gaussian distribution.
        minimum = the minimum value of the probability density to be created.
        maximum = the maximum value of the probability density to be created.
        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = np.linspace(minimum, maximum, num=100, endpoint=True)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf = (np.exp(-(fp - mu) ** 2 / (2 * sigma ** 2))).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return_array = np.interp(interpolator, xp, fp)

        return return_array

    def sin_gauss_range(self, mu, sigma, minimum, maximum):
        """
        Obj method used to create an array corresponding to a probability density of a gaussian distribution multiplied
        by a sine function. (i.e. y=sin(x)*exp(-(x-mu)^2/(2*sigma^2))

        Parameters
        ----------
        mu = The mean of the gaussian distribution.
        sigma = The standard deviation of the gaussian distribution.
        minimum = the minimum value of the probability density to be created.
        maximum = the maximum value of the probability density to be created.
        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = np.linspace(minimum, maximum, num=100, endpoint=True)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf = (np.sin(fp) * np.exp(-(fp - mu) ** 2 / (2 * sigma ** 2))).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return_array = np.interp(interpolator, xp, fp)

        return return_array

    def linear_peak(self, minimum, mid, maximum):
        """
        Obj method used to create an array corresponding to a probability density that rises linearly from zero at a
        given minimum bound to a maximum probability at a given mid value, before dropping linearly back to zero at a
        given maximum bound. (i.e. for values less than the mid value, (x-minimum)/(mid-minimum), and for values
        greater than the mid value, (maximum-x)/(maximum-mid).)

        Parameters
        ----------
        minimum = The minimum value of the distribution where probability starts at zero before rising linearly towards
            the mid value.
        mid = The value of maximum probability between the values of minimum and maximum.
        maximum = The maximum value of the distribution where probability ends at zero, dropping linearly from
            the mid value.
        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = np.linspace(minimum, maximum, num=100, endpoint=True)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf1 = (1. / (mid - minimum) * fp[fp < mid] - minimum / (mid - minimum))
        cdf2 = (-1. / (maximum - mid) * fp[fp >= mid] + maximum / (maximum - mid))
        # Combines the two non-continuous halves of our distribution.
        cdf = np.concatenate([cdf1, cdf2]).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return_array = np.interp(interpolator, xp, fp)

        return return_array

    def power_knee_divot(self, alpha_bright, h_max, h_min=1., h_break=None, alpha_faint=None, contrast_ratio=1.):
        """
        Obj method used to create an array of H-magnitudes from a so-called single power-law, knee,
        or divot H-magnitude distribution.

        Parameters
        ----------
        alpha_bright = The slope of the bright section of the power law.
        h_max = Maximum H value generated.
        h_min = (optional) Minimum H value generated (default=1).
        h_break = (optional) The H-magnitude where the function transitions from the bright distribution to the
            faint distribution (default=None).
        alpha_faint = (optional) The slope of the faint section of the power law (default=None).
        contrast_ratio = (optional) The contrast is the ratio of the number of bright objects at h_break to the
            number of faint objects at h_break. Contrast defines the "size" of the divot.
            e.g. A contrast of 1  means there is no divot, while a contrast of 10 means
            there are 10 times as many objects in the bright population as there are in the faint population for an
            H-magnitude equal to h_break. (default=1).

        When provided a slope alpha and a faint-side maximum H-magnitude (h_max), an H-magnitude is drawn randomly
        from the distribution:
            dN/dH = k*ln(10)*alpha*10^(alpha*H)
        in the range h_min = 1 to h_max. Specify h_min to change the bright-end.

        Specifying an h_break and alpha_faint will draw from a knee distribution.

        Specifying an h_break, alpha_faint and contrast_ratio will draw from a divot distribution as in
        Shankman et al. 2013

        e.g.

        ---Single Power Law---
        draw_h(0.8,13)
        will draw an H-magnitude from the appropriate distribution such that H [1,13]
        draw_h(0.8,13,h_min=5)
        will draw an H-magnitude such that H [5,13]

        ---Knee---
        To draw from a knee distribution specify hbreak and alpha_faint
        draw_h(0.8, 13, hbreak=9, alpha_faint = 0.5)
        This will draw an H-magnitude from a distribution that breaks at H=9 from a slope of 0.8 to a slope of 0.5.
        hmin can also be specified here.

        ---Divot---
        To draw from a divot (see Shankman et al 2013), specify hbreak, alpha_faint, and the contrast value.
        Contrasts should be > 1. hmin can also be specified.
        draw_h(0.8, 13, hbreak=9, alpha_faint = 0.5, contrast = 23)
        """

        # Simple error handling for silly input
        if h_max < h_break:
            raise ValueError('h_max must be greater than h_break')
        if h_max < h_min:
            raise ValueError('h_max must be greater than h_min')
        if h_break is not None and h_break < h_min:
            raise ValueError('h_break must be greater than h_min')

        # Avoid singularity for alpha_bright = 0
        alpha_bright = 0.0000000001 if alpha_bright == 0 else alpha_bright

        # Set alpha_faint to alpha_bright for the case of a single power-law
        alpha_faint = alpha_bright if alpha_faint is None else alpha_faint

        # Avoid singularity for alpha_faint = 0
        alpha_faint = 0.0000000001 if alpha_faint == 0 else alpha_faint

        # Set h_break to be the maximum H for the case of a single power-law
        h_break = h_max if h_break is None else h_break

        # Generate the bright side of the size distribution
        fp_b = np.linspace(h_min, h_break, num=100, endpoint=True)
        xp_b = 10 ** (alpha_bright * fp_b) * alpha_bright * np.log(10)

        if h_break is not None:
            # Generate the faint side of the size distribution
            fp_f = np.linspace(h_break, h_max, num=100, endpoint=True)
            xp_f = 10 ** (alpha_faint * fp_f) * alpha_faint * np.log(10)
        else:
            fp_f = []
            xp_f = []

        # Coefficient to merge the faint and bright xp arrays
        xp_f = xp_f * ((alpha_bright * 10 ** (alpha_bright * h_break)) / (alpha_faint * 10 ** (alpha_faint * h_break))
                       / contrast_ratio)

        # Merge the x axes into one array
        fp = np.concatenate([fp_b, fp_f])

        # Merge them and then accumulate and normalize
        cdf = np.concatenate([xp_b, xp_f]).cumsum()
        xp = cdf / max(cdf)

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return_array = np.interp(interpolator, xp, fp)

        return return_array

    def increment_index(self):
        """
        Used to increment the index number after the current array value has been used.

        This method also re-draws the Obj instances when the index exceeds the number of elements in the incremented
        instance.
        """

        self.index = self.index + 1

        if self.index > self.size:
            initialize()

    def return_values(self):

        a = self.a[self.index]
        e = self.e[self.index]
        i = self.i[self.index]
        node = self.node[self.index]
        peri = self.peri[self.index]
        m = self.m[self.index]
        epoch = self.epoch[self.index]
        h = self.h[self.index]
        color = self.color[self.index, :]
        gb = self.gb[self.index]
        ph = self.ph[self.index]
        period = self.period[self.index]
        amp = self.amp[self.index]
        resamp = self.resamp[self.index]
        comments = self.comments

        return a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments


def set_seed(seed):
    """
    This function initializes the numpy random number generator using a given seed.
    It assigns this generator to a global variable so it is accessible in other parts of this module.
    All subsequent calls for a random number must be made from this generator instance.
        e.g. rnd.uniform(0,1,10**6)

    Calling np.random.uniform(0,1,10**6) directly operates off of a different generator which uses a different
    seed, not set by the user.

    Parameters
    ----------
    seed = User chosen seed value for the random number generator.
    """
    global rnd_gen

    rnd_gen = np.random.RandomState(seed)  # Creates generator instance.


def initialize():
    """
    Initializes Obj class objects containing the arrays of a given size before filling them with the population
    distributions given in a "draw_" method.
    """

    size = 10 ** 6

    # Initialize class instances and fill their arrays.
    global obj1, obj2, obj3
    obj1 = Obj(size)
    obj1.draw_symm()  # Symmestric resonance

    obj2 = Obj(size)
    obj2.draw_asymm_lead()  # Asymmetric leading resonance

    obj3 = Obj(size)
    obj3.draw_asymm_trail()  # Asymmetric trailing resonance

    # import pylab as plt
    # fig = plt.figure(figsize=(6, 6))
    # plt.plot(Var.xp_resamp, Var.fp_resamp, '.k', markersize=1)
    # plt.show()
    # plt.plot(np.append(cdf1, cdf2))
    # plt.show()
    # plt.plot(Var.xp_resamp)
    # plt.show()
    # plt.plot(Var.resamp, '.k', markersize=1)
    # plt.show()


def gimeobj():
    """
    gimeobj is the function that determines which object instance to draw values from, and then returns those values
    to the Driver module.

    Output
    ------
    a - semimajor axis
    e - eccentricity
    i - inclination [rad]
    node - longitude of the ascending node [rad]
    peri - argument of pericenter [rad]
    m - mean anomaly (M) [rad]
    h - H magnitude in field of survey
    color - array of colour conversion values
    gb - opposition surge effect
    ph - phase angle
    period - light curve period
    amp - light curve amplitude
    resamp - resonant amplitude
    comments - comment field to indicate source population
    """

    if rnd_gen.uniform() < 0.3:  # Symmetric
        # Assign values from the object instance to be returned to Driver, and then increment the array index.
        a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments = obj1.return_values()
        obj1.increment_index()
    else:
        if rnd_gen.uniform() < 0.5:  # Asymmetric leading
            # Assign values from the object instance to be returned to Driver, and then increment the array index.
            a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments = obj2.return_values()
            obj2.increment_index()
        else:  # Asymmetric trailing
            # Assign values from the object instance to be returned to Driver, and then increment the array index.
            a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments = obj3.return_values()
            obj3.increment_index()

    return a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments

    # The call to SurveySubs.detos1 for reference
    # SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)
