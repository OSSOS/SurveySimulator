"""GiMeObj is a module used by Driver to run the Survey Simulator. This module contains all of the required functions
and methods to generate and return a series of objects to Driver which then feeds them into the fortran code
comprising the computationally intensive parts of the Survey Simulator.

GiMeObj is intended to be modified by the user to create distributions of objects representing theoretical
populations of objects in the distant solar system. These theoretical objects are then fed into Driver to evaluate if
they would have been detected by surveys detailed in the appropriate folder of the survey simulator."""
import numpy as np
import matplotlib.pyplot as plt  # Only needed if using the diagnostic plotting methods.

rnd_gen = None  # Creates empty variable to hold the random number generator once the seed is set.

# Creates empty global variables to house the object class instances generated in initialize().
obj1 = None
obj2 = None
obj3 = None
obj4 = None
obj5 = None


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

    def draw_distribution(self):
        """
        Initializes an Obj class object used to contain the arrays of a given size filled with the population
        distributions of the various orbital parameters for a single orbital population (e.g. the plutinos).

        This distribution represents a more complicated example than provided in the "Ring" module. While the
        distributions presented here are not intended to accurately reflect real-world TNO populations, they serve to
        illustrate the various methods that can be used to generate parameter distributions.
        """

        self.comments = "Distribution"

        self.a = self.uniform(47.8, 48.8)  # units in AU
        self.e = self.gauss_range(0.18, 0.06, 0, 1)
        self.i = self.sin_gauss_range(np.deg2rad(0), np.deg2rad(16), np.deg2rad(0), np.deg2rad(90))  # units in radians

        self.node = self.uniform(0, 2*np.pi)  # random value between 0 and 2pi # units in radians
        self.m = self.uniform(0, 2*np.pi)   # random value between 0 and 2pi # units in radians
        self.epoch = self.constant(2456293.5)

        # Import constants characterizing the H-magnitude distribution to be generated.
        alpha_f, alpha_b, hmax, hmin, hbreak, contrast = \
            np.genfromtxt('H_dist_vars.txt', unpack=True, comments='#')
        self.h = self.power_knee_divot(alpha_b, hmax, hmin, hbreak, alpha_f, contrast)

        self.gb = self.constant(-0.12)  # Opposition surge effect
        self.ph = self.constant(0.00)  # Initial phase of lightcurve  # radians
        self.period = self.constant(0.60)  # Period of lightcurve  # units in days
        self.amp = self.constant(0.00)  # Amplitude of lightcurve (peak-to-peak)

        self.resamp = self.linear_peak(np.deg2rad(0), np.deg2rad(5), np.deg2rad(10))  # units in radians

        # choose value of phi sinusoidally from resamp
        phi = np.pi + np.sin(self.uniform(0, 2*np.pi))*self.resamp
        lambda_n = np.deg2rad(333.5078)  # Neptune's mean longitude on 1 Jan 2013

        self.peri = (1 / 1 * (phi - 2 * self.m) - self.node + lambda_n) % (2 * np.pi)  # units in radians

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

    def read_file(self, f_name, column):
        """
        Obj method used to read in data from a file and assign it in an array to one of the Obj variables.

        read_file must be the first methods called to assign arrays in the draw_X method where it is used.

        Once data has been read in from a file this method duplicates this set of data until it has reached the size
        given at object instantiation, or greater. It then changes the instance variable "size" to match the length
        of the generated array. As such all methods called after this point will have the same length.

        If multiple calls are made to read_file then all input data must have the same number of entries for each
        call. If not, the array indexing will not function as expected, and some variables will not have
        corresponding values for the full data sets.

        If a data set larger than that provided by the input file is required, it is left to the user to decide how
        to generate this expanded data set based on their individual needs. Looping over the same repeated values
        read in from a file can yield unexpected results from the survey simulator.

        Parameters
        ----------
        f_name = A string giving the name of the input file.
        column = The column number where the data of interest can be found, counting from 0.
        """
        file_array = np.genfromtxt(f_name, usecols=column, unpack=True, comments='#')

        return_array = file_array

        # Appends the set of data returned from the file onto the end of the return array until the return array
        # equals or exceeds the instance "size" variable.
        while return_array.size < self.size:
            return_array = np.concatenate([return_array, file_array])

        self.size = return_array.size  # Sets the instance "size" variable to this new maximum length.

        return return_array

    def blur(self, attribute, fraction):
        """
        Method to increase or decrease values in a given variable up to a given fraction of that variable. This
        "blurring" alters all values in an array by a uniformly generated random value.

        This method is intended to create "extra" values from a limited data set read in from a file. it will help to
        avoid unintended statistical effects for repeated values read into the survey simulator. This specific
        implementation is not necessarily appropriate for all data sets read in, but instead serves as a simple
        example for the types of operations that would need to be performed. The exact form that this "blurring" should
        take is left to be developed by the user for their specific needs.

        It is recommended that the diagnostic .plot() method be used to examine any generated distributions to ensure
        the new arrays behave as expected.

        Parameters
        ----------
        variable - The variable to be "blurred"". Given as a string.
        fraction - The fraction used to alter the variable. Fraction must be between 0 and 1, inclusive.
        """

        if fraction < 0 or fraction > 1:  # Error checking for out of bounds fraction choice.
            raise ValueError("Fraction must be between 0 and 1, inclusive.")

        source = getattr(self, attribute)  # Store chosen array as variable "source" for ease of use.
        change = source*fraction  # Maximum increase or decrease.

        # Generate new "blurred" array by generating an array (of the same size as the source array)
        # of uniform random values on the interval [-change, change) before adding this new array element-wise
        # to the original source array.
        blurred = rnd_gen.uniform(-change, change, source.size) + source

        # Set original attribute to the new array generated from the source attribute.
        setattr(self, attribute, blurred)

    def blur_e(self, attribute, fraction):
        """
        A special case of the blur method defined above. See that method for full documentation.

        Eccentricity needs special consideration because the above "blurring" method can easily yield values over 1,
        which correspond to orbits that are unlikely to be of interest.

        Handling results where eccentricity is greater than 1 is extremely complicated. In lieu of a complicated and
        computationally intensive procedure, that is unlikely to be of assistance to most users, this method instead
        provides a simple way to handle these problematic values. Eccentricity values where e*(1+fraction) >= 1 are
        not "blurred" and instead set to their original value. This results in not ideal behavior for high e values,
        but is a functional solution for this simple example method.

        The user will want to handle high e values with more care when creating their own methods.

        Parameters
        ----------
        variable - The variable to be "blurred"". Given as a string.
        fraction - The fraction used to alter the variable. Fraction must be between 0 and 1, inclusive.
        """

        if fraction < 0 or fraction > 1:  # Error checking for out of bounds fraction choice.
            raise ValueError("Fraction must be between 0 and 1, inclusive.")

        source = getattr(self, attribute)  # Store chosen array as variable "source" for ease of use.
        change = source * fraction  # Maximum increase or decrease.

        # Generate new "blurred" array by generating an array (of the same size as the source array)
        # of uniform random values on the interval [-change, change) before adding this new array element-wise
        # to the original source array.
        blurred = rnd_gen.uniform(-change, change, source.size) + source

        # Checks for high e values that will be out of the range of interest (equal to or greater than 1)
        # and sets them to their original value.
        for index in range(0, source.size):
            if source[index]+change[index] >= 1:
                blurred[index] = source[index]

        # Set original attribute to the new array generated from the source attribute.
        setattr(self, attribute, blurred)

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
        Obj method that creates an array of the appropriate size filled with random variables on the given interval.
        (default interval [0,1).)

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

        This method has the added restriction of placing a maximum and minimum value on the generated distribution.

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

        This method has the added restriction of placing a maximum and minimum value on the generated distribution.

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

        This method also re-draws all of the Obj instances when the index exceeds the number of elements in the
        incremented instance.

        When using numerous Obj instances it would be of use to check and only re-initialize specific instances that
        have used up all the values therein as opposed to re-initializing all instances as this method does now. Such
        an operation adds program complexity, but for large numbers of different instances will save computation
        time. The trade-off between complexity and optimization is left to the user.
        """

        self.index += 1

        if self.index > self.size - 1:
            initialize()

    def return_values(self):
        """
        Method used to return the values from the arrays at the current index to the main gimeobj function. This
        method selects the values to be returned to Driver for evaluation in the fortran survey simulator code.
        """

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

    # Methods below this point are intended for diagnostic purposes, not for generating populations or the basic
    # operations of the module overall.
    def plot(self, *variables):
        """
        Method that will plot GiMeObj output variables stored in the specified Obj instance.
        The parameter arrays are plotted as the parameter value vs position in the array (index).

        Accepts a variable number of arguments. Arguments must be supplied as strings.
        (e.g. instance.plot("a", "e") would plot just the semi-major axis and eccentricity arrays.)
        No arguments defaults to plotting all variables.

        Once all plots have been generated and displayed this method terminates the program, never returning to Driver.

        Intended for diagnostic purposes.

        Parameter Options
        -----------------
        "a" - Semi-major Axis [AU]
        "e" - Eccentricity
        "i" or "inc" - Inclination [°]
        "node" - Longitude of the Ascending Node (Ω) [°]
        "peri" - Argument of Pericenter (ω) [°]
        "m" or "M" - Mean Anomaly (M) [°]
        "epoch" - Epoch [JD]
        "h" or "H" - H-Magnitude
        "color" - 9 graphs representing the color information in the 2D color array.
        "gb" - Opposition Surge Effect
        "ph" - Phase Angle [°]
        "period" - Light Curve Period [days]
        "amp" - Light Curve Amplitude
        "resamp" - Resonant Amplitude [°]\
        "all" or "All" - Plots all variables.
        """
        flag = True  # Flag used to determine if no arguments were passed to the plot method.

        for select in variables:
            flag = False

            if select == "a":
                self.plot_a()
            elif select == "e":
                self.plot_e()
            elif select == "i" or select == "inc":
                self.plot_i()
            elif select == "node":
                self.plot_node()
            elif select == "peri":
                self.plot_peri()
            elif select == "m" or select == "M":
                self.plot_m()
            elif select == "epoch":
                self.plot_epoch()
            elif select == "h" or select == "H":
                self.plot_h()
            elif select == "color":
                self.plot_color()
            elif select == "gb":
                self.plot_gb()
            elif select == "ph":
                self.plot_ph()
            elif select == "period":
                self.plot_period()
            elif select == "amp":
                self.plot_amp()
            elif select == "resamp":
                self.plot_resamp()
            elif select == "all" or select == "All":
                self.plot_all()
            else:
                raise ValueError("Invalid variable argument.")

        if flag:
            self.plot_all()

        exit()  # Terminates execution after all plots have been generated.

    def plot_a(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.a, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Semi-major Axis [AU]")
        plt.show()

    def plot_e(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.e, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Eccentricity")
        plt.show()

    def plot_i(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.i), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Inclination [°]")
        plt.show()

    def plot_node(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.node), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Longitude of the Ascending Node (Ω) [°]")
        plt.show()

    def plot_peri(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.peri), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Argument of Pericenter (ω) [°]")
        plt.show()

    def plot_m(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.m), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Mean Anomaly (M) [°]")
        plt.show()

    def plot_epoch(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.epoch, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Epoch [JD]")
        plt.show()

    def plot_h(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.h, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("H-Magnitude")
        plt.show()

    def plot_color(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 0], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, g-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 1], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, r-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 2], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, i-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 3], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, z-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 4], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, u-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 5], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, V-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 6], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, B-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 7], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, R-x")
        plt.show()

        plt.figure(figsize=(15, 12))
        plt.plot(self.color[:, 8], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Color, I-x")
        plt.show()

        # plt.figure(figsize=(15, 12))
        # plt.plot(self.color[:, 9], color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        # plt.ylabel("Color, extra")
        # plt.show()

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

    def plot_gb(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.gb, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Opposition Surge Effect")
        plt.show()

    def plot_ph(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.ph), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Phase Angle [°]")
        plt.show()

    def plot_period(self):
        plt.figure(figsize=(15, 12))
        plt.plot(self.period, color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Light Curve Period [days]")
        plt.show()

    def plot_amp(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.amp), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Light Curve Amplitude")
        plt.show()

    def plot_resamp(self):
        plt.figure(figsize=(15, 12))
        plt.plot(np.rad2deg(self.resamp), color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
        plt.ylabel("Resonant Amplitude [°]")
        plt.show()

    def plot_all(self):
        self.plot_a()
        self.plot_e()
        self.plot_i()
        self.plot_node()
        self.plot_peri()
        self.plot_m()
        self.plot_epoch()
        self.plot_h()
        self.plot_color()
        self.plot_gb()
        self.plot_ph()
        self.plot_period()
        self.plot_amp()
        self.plot_resamp()


def set_seed(seed):
    """
    This function initializes the numpy random number generator using a given seed.
    It assigns this generator to a global variable so it is accessible in other parts of this module.
    All subsequent calls for a random number must be made from this generator instance.
        e.g. rnd_gen.uniform(0,1,10**6)

    Calling np.random.uniform(0,1,10**6) directly operates off of a different generator instance which uses a different
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
    global obj1
    obj1 = Obj(size)
    obj1.draw_distribution()  # Generate distribution.

    # Plotting method used to examine values contained in object arrays.
    # obj1.plot("all")

def gimeobj():
    """
    gimeobj is the function that determines which object instance to draw values from, and then returns those values
    to the Driver module.

    Output
    ------
    a - semi-major axis [AU]
    e - eccentricity
    i - inclination [rad]
    node - longitude of the ascending node (Ω) [rad]
    peri - argument of pericenter (ω) [rad]
    m - mean anomaly (M) [rad]
    epoch - Time of elements [JD]
    h - H magnitude in field of survey
    color - array of colour conversion values
    gb - opposition surge factor G, Bowell formalism
    ph - phase of light curve at epoch jday [rad]
    period - light curve period [days]
    amp - light curve amplitude [mag]
    resamp - resonant amplitude [rad]
    comments - comment field to indicate source population
    """

    # Assign values from the object instance to be returned to Driver.
    a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments = obj1.return_values()
    # Increment the array index value and check if the arrays need to be redrawn.
    obj1.increment_index()

    # Return values to Driver.
    return a, e, i, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp, comments

    # The call to SurveySubs.detos1 for reference
    # SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)
