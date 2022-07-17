"""
A module of methods that useful for creating distributions of values used to model KBO orbital distributions.
"""
import numpy


class Distributions(object):
    """
    The Distribution class holds a number of methods used to generate orbit parameter distributions appropriate for use as input
    models to the Survey Simulator.
    """

    def __init__(self, seed=None, size=10**6):
        """
        Args:
            seed (None or int): the seed used to intialize the random number generator.
            size (int): how many samples to generate with each pass.

        The advantage of Distribution is that a sample of 'size' is created at each call and then read from for each object generated,
        this is faster than calling each distribution function for every object needed.
        """
        self.seed = seed
        self.size = size
        self.rnd_gen = numpy.random.default_rng(self.seed)

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
        change = source * fraction  # Maximum increase or decrease.

        # Generate new "blurred" array by generating an array (of the same size as the source array)
        # of uniform random values on the interval [-change, change) before adding this new array element-wise
        # to the original source array.
        blurred = self.rnd_gen.uniform(-change, change, self.size) + source

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
        blurred = self.rnd_gen.uniform(-change, change, self.size) + source

        # Checks for high e values that will be out of the range of interest (equal to or greater than 1)
        # and sets them to their original value.
        for index in range(0, self.size):
            if source[index] + change[index] >= 1:
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

        return_array = numpy.full(self.size, constant)

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

        return_array = self.rnd_gen.uniform(minimum, maximum, self.size)

        return return_array

    def normal(self, mu, sigma):
        """
        Obj method that creates an array of the appropriate size filled with random variables in a gaussian
        distribution.

        Parameters
        ----------
        mu = The mean of the gaussian distribution.
        sigma = The standard deviation of the gaussian distribution.
        """

        return_array = self.rnd_gen.normal(mu, sigma, self.size)

        return return_array

    def truncated_normal(self, mu, sigma, minimum, maximum):
        """
        Probability density of a gaussian distribution. (i.e. y=exp(-(x-mu)^2/(2*sigma^2))

        This method has the added restriction of placing a maximum and minimum value on the generated distribution.

        Args:
            mu (float): The mean of the gaussian distribution.
            sigma (float): The standard deviation of the gaussian distribution.
            minimum (float): the minimum value of the probability density to be created.
            maximum (float): the maximum value of the probability density to be created.
        Returns:
            numpy.array: random distribution of values following a truncated_normal distribution.
        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = numpy.linspace(minimum, maximum, num=100)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf = (numpy.exp(-(fp - mu) ** 2 / (2 * sigma ** 2))).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = self.rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return numpy.interp(interpolator, xp, fp)

    def truncated_sin_normal(self, mu, sigma, minimum, maximum):
        """
        Obj method used to create an array corresponding to a probability density of a gaussian distribution multiplied
        by a sine function. (i.e. y=sin(x)*exp(-(x-mu)^2/(2*sigma^2)) (follows from Brown 2001)

        This method has the added restriction of placing a maximum and minimum value on the generated distribution.

        Args:
            mu (float): The mean of the gaussian distribution, in radians
            sigma (float): The standard deviation of the gaussian distribution, in radians.
            minimum (float): the minimum value of the probability density to be created, in radians.
            maximum (flaot): the maximum value of the probability density to be created, in radians.
        Returns:
            numpy.array: resulting random samples.
        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = numpy.linspace(minimum, maximum, num=100)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf = (numpy.sin(fp) * numpy.exp(-(fp - mu) ** 2 / (2 * sigma ** 2))).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = self.rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return numpy.interp(interpolator, xp, fp)

    def linear_peak(self, minimum, mid, maximum):
        """
        Generate an array corresponding to a probability density that rises linearly from zero at a
        given minimum bound to a maximum probability at a given mid value, before dropping linearly back to zero at a
        given maximum bound. (i.e. for values less than the mid value, (x-minimum)/(mid-minimum), and for values
        greater than the mid value, (maximum-x)/(maximum-mid).)

        Args:
            minimum (flaot): minimum value of the distribution where probability starts at zero before rising linearly
            towards the mid value.
            mid (float): value of maximum probability between the values of minimum and maximum.
            maximum (float): maximum value of the distribution where probability ends at zero, dropping linearly
            from the mid value.
        Returns:
            numpy.array: resulting random samples.

        """

        # Create an array, of size num, evenly spaced values from minimum to maximum including the maximum value.
        fp = numpy.linspace(minimum, maximum, num=100)

        # The desired probability density is evaluated over the desired range of values given in the array fp.
        # This creates a new array which is summed over using the cumulative sum (cumsum()) function.
        # This array functions as our cumulative distribution function (cdf).
        cdf1 = (1. / (mid - minimum) * fp[fp < mid] - minimum / (mid - minimum))
        cdf2 = (-1. / (maximum - mid) * fp[fp >= mid] + maximum / (maximum - mid))
        # Combines the two non-continuous halves of our distribution.
        cdf = numpy.concatenate([cdf1, cdf2]).cumsum()
        xp = cdf / max(cdf)  # Normalizes the cdf

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = self.rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return numpy.interp(interpolator, xp, fp)

    def power_knee_divot(self, alpha_bright, h_max, h_min=1., h_break=None, alpha_faint=None, contrast_ratio=1.):
        """
        Obj method used to create an array of H-magnitudes from a so-called single power-law, knee,
        or divot H-magnitude distribution.

        Args:
        alpha_bright (float): slope of the bright section of the power law.
        h_max (float): maximum H value generated.
        h_min (float): minimum H value generated (default=1).
        h_break (float): H-magnitude where the function transitions from the bright distribution to the
        faint distribution (default=None, no break).
        alpha_faint (float): slope of the faint section of the power law (default=None).
        contrast_ratio (float): The contrast is the ratio of the number of bright objects at h_break to the
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
        To draw from a knee distribution specify h_break and alpha_faint
        draw_h(0.8, 13, h_break=9, alpha_faint = 0.5)
        This will draw an H-magnitude from a distribution that breaks at H=9 from a slope of 0.8 to a slope of 0.5.
        h_min can also be specified here.

        ---Divot---
        To draw from a divot (see Shankman et al 2013), specify h_break, alpha_faint, and the contrast value.
        Contrasts should be > 1. h_min can also be specified.
        draw_h(0.8, 13, h_break=9, alpha_faint = 0.5, contrast = 23)
        """

        # Simple error handling for silly input
        if h_max < h_break:
            raise ValueError('h_max must be greater than h_break')
        if h_max < h_min:
            raise ValueError('h_max must be greater than h_min')
        if h_break is not None and h_break < h_min:
            raise ValueError('h_break must be greater than h_min')

        # Avoid singularity for alpha_bright = 0
        alpha_bright = 1e-10 if alpha_bright == 0 else alpha_bright

        # Set alpha_faint to alpha_bright for the case of a single power-law
        alpha_faint = alpha_bright if alpha_faint is None else alpha_faint

        # Avoid singularity for alpha_faint = 0
        alpha_faint = 1e-10 if alpha_faint == 0 else alpha_faint

        # Set h_break to be the maximum H for the case of a single power-law
        h_break = h_max if h_break is None else h_break

        # Generate the bright side of the size distribution
        fp_b = numpy.linspace(h_min, h_break, num=100)
        xp_b = 10 ** (alpha_bright * fp_b) * alpha_bright * numpy.log(10)

        if h_break is not None:
            # Generate the faint side of the size distribution
            fp_f = numpy.linspace(h_break, h_max, num=100)
            xp_f = 10 ** (alpha_faint * fp_f) * alpha_faint * numpy.log(10)
        else:
            fp_f = []
            xp_f = []

        # Coefficient to merge the faint and bright xp arrays
        xp_f = xp_f * ((alpha_bright * 10 ** (alpha_bright * h_break)) / (alpha_faint * 10 ** (alpha_faint * h_break))
                       / contrast_ratio)

        # Merge the x axes into one array
        fp = numpy.concatenate([fp_b, fp_f])

        # Merge them and then accumulate and normalize
        cdf = numpy.concatenate([xp_b, xp_f]).cumsum()
        xp = cdf / max(cdf)

        # Generates an array (of a size given in our instance initialization) of uniformly distributed random values
        # ranging from the minimum value of of our normalized cdf (xp[0] which is >=0) to the maximum value (1).
        interpolator = self.rnd_gen.uniform(xp[0], 1, self.size)

        # Calls the interpolation function to generate a randomized array, corresponding to the probability
        # density, by interpolating between the given points.
        return numpy.interp(interpolator, xp, fp)
