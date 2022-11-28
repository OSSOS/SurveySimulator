from copy import deepcopy
import warnings

from astropy.units import Quantity

__author__ = "David Rusk <drusk@uvic.ca>"

import math
from astropy import wcs as astropy_wcs
from astropy import units
import numpy
import logging


PI180 = 57.2957795130823208767981548141052


class WCS(astropy_wcs.WCS):
    def __init__(self, header):
        """
        Create the bits needed for working with sky2xy
        """
        astropy_header = deepcopy(header)
        del (astropy_header['PV*'])
        with warnings.catch_warnings(record=True):
            warnings.resetwarnings()
            warnings.simplefilter(
                "ignore", astropy_wcs.FITSFixedWarning, append=True)
            super(WCS, self).__init__(astropy_header)
        self.header = header

    @property
    def cd(self):
        """
        CD Rotation matrix values.
        """
        return parse_cd(self.header)

    @property
    def dc(self):
        """
        CD Rotation matrix INVERTED i.e.  []^-1
        """
        return numpy.linalg.inv(self.cd)
        # return numpy.array(numpy.mat(self.cd).I)

    @property
    def pv(self):
        """
        Array of PV keywords used for hi-odered astrogwyn mapping
        """
        return parse_pv(self.header)

    @property
    def crpix1(self):
        """
        1st reference coordinate
        """
        return self.header['CRPIX1']

    @property
    def crpix2(self):
        """
        2nd reference coordinate
        """
        return self.header['CRPIX2']

    @property
    def crval1(self):
        """
        Reference Coordinate of 1st reference pixel
        """
        return self.header['CRVAL1']

    @property
    def crval2(self):
        """
        Reference Coordinate of 2nd reference pixel
        """
        return self.header['CRVAL2']

    @property
    def nord(self):
        """
        The order of the PV fit, provided by astgwyn
        """
        return self.header['NORDFIT']

    def xy2sky(self, x, y, usepv=True):
        if usepv:
            try:
                return xy2skypv(x=numpy.array(x), y=numpy.array(y),
                                crpix1=self.crpix1,
                                crpix2=self.crpix2,
                                crval1=self.crval1,
                                crval2=self.crval2,
                                cd=self.cd,
                                pv=self.pv,
                                nord=self.nord)
            except Exception as ex:
                logging.warning("Error {} {}".format(type(ex), ex))
                logging.warning("Reverted to CD-Matrix WCS.")
        xy = numpy.array([x, y]).transpose()
        pos = self.wcs_pix2world(xy, 1).transpose()
        return pos[0] * units.degree, pos[1] * units.degree

    def sky2xy(self, ra, dec, usepv=True):
        if isinstance(ra, Quantity):
            ra = ra.to(units.degree).value
        if isinstance(dec, Quantity):
            dec = dec.to(units.degree).value
        try:
            if usepv:
                return sky2xypv(ra=ra,
                                dec=dec,
                                crpix1=self.crpix1,
                                crpix2=self.crpix2,
                                crval1=self.crval1,
                                crval2=self.crval2,
                                dc=self.dc,
                                pv=self.pv,
                                nord=self.nord)
        except Exception as ex:
            logging.warning("sky2xy raised exception: {0}".format(ex))
            logging.warning("Reverted to CD-Matrix WCS to convert: {0} {1} ".format(ra, dec))
        pos = self.wcs_world2pix([[ra, dec], ], 1)
        return pos[0][0], pos[0][1]


def sky2xypv(ra, dec, crpix1, crpix2, crval1, crval2, dc, pv, nord, maxiter=300):
    """
    Transforms from celestial coordinates to pixel coordinates to taking
    non-linear distortion into account with the World Coordinate System
    FITS keywords as used in MegaPipe.

    For the inverse operation see xy2sky.

    Reference material:
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/megapipe/docs/CD_PV_keywords.pdf

    Args:
      ra: float
        Right ascension
      dec: float
        Declination
      crpix1: float
        Tangent point x, pixels
      crpix2: float
        Tangent point y, pixels
      crval1: float
        Tangent point RA, degrees
      crval2:
        Tangent point Dec, degrees
      dc: 2d array
        Expresses the scale, the rotation and any possible skew of the image
        with respect to the sky.  It is the inverse of the cd matrix in xy2sky.
      pv: 2d array
      nord: int
        order of the fit

    Returns:
      x, y: float
        Pixel coordinates
    """
    if math.fabs(ra - crval1) > 100:
        if crval1 < 180:
            ra -= 360
        else:
            ra += 360

    ra /= PI180
    dec /= PI180

    tdec = math.tan(dec)
    ra0 = crval1 / PI180
    dec0 = crval2 / PI180
    ctan = math.tan(dec0)
    ccos = math.cos(dec0)

    traoff = math.tan(ra - ra0)
    craoff = math.cos(ra - ra0)
    etar = (1 - ctan * craoff / tdec) / (ctan + craoff / tdec)
    xir = traoff * ccos * (1 - etar * ctan)
    xi = xir * PI180
    eta = etar * PI180

    if nord < 0:
        # The simple solution
        x = xi
        y = eta
    else:
        # Reverse by Newton's method
        tolerance = 0.001 / 3600

        # Initial guess
        x = xi
        y = eta
        iteration = 0

        converged = False
        while not converged:
            assert nord >= 0

            # Estimates
            f = pv[0][0]
            g = pv[1][0]

            # Derivatives
            fx = 0
            fy = 0
            gx = 0
            gy = 0

            if nord >= 1:
                r = math.sqrt(x ** 2 + y ** 2)
                f += pv[0][1] * x + pv[0][2] * y + pv[0][3] * r
                g += pv[1][1] * y + pv[1][2] * x + pv[1][3] * r
                fx += pv[0][1] + pv[0][3] * x / r
                fy += pv[0][2] + pv[0][3] * y / r
                gx += pv[1][2] + pv[1][3] * x / r
                gy += pv[1][1] + pv[1][3] * y / r

                if nord >= 2:
                    x2 = x ** 2
                    xy = x * y
                    y2 = y ** 2

                    f += pv[0][4] * x2 + pv[0][5] * xy + pv[0][6] * y2
                    g += pv[1][4] * y2 + pv[1][5] * xy + pv[1][6] * x2
                    fx += pv[0][4] * 2 * x + pv[0][5] * y
                    fy += pv[0][5] * x + pv[0][6] * 2 * y
                    gx += pv[1][5] * y + pv[1][6] * 2 * x
                    gy += pv[1][4] * 2 * y + pv[1][5] * x

                    if nord >= 3:
                        x3 = x ** 3
                        x2y = x2 * y
                        xy2 = x * y2
                        y3 = y ** 3

                        f += pv[0][7] * x3 + pv[0][8] * x2y + pv[0][9] * xy2 + pv[0][10] * y3
                        g += pv[1][7] * y3 + pv[1][8] * xy2 + pv[1][9] * x2y + pv[1][10] * x3
                        fx += pv[0][7] * 3 * x2 + pv[0][8] * 2 * xy + pv[0][9] * y2
                        fy += pv[0][8] * x2 + pv[0][9] * 2 * xy + pv[0][10] * 3 * y2
                        gx += pv[0][8] * y2 + pv[1][9] * 2 * xy + pv[1][10] * 3 * x2
                        gy += pv[1][7] * 3 * y2 + pv[0][8] * 2 * xy + pv[1][9] * x2

            f -= xi
            g -= eta
            dx = (-f * gy + g * fy) / (fx * gy - fy * gx)
            dy = (-g * fx + f * gx) / (fx * gy - fy * gx)
            x += dx
            y += dy

            if math.fabs(dx) < tolerance and math.fabs(dy) < tolerance:
                converged = True

            iteration += 1

            if iteration > maxiter:
                break

    xp = dc[0][0] * x + dc[0][1] * y
    yp = dc[1][0] * x + dc[1][1] * y

    x = xp + crpix1
    y = yp + crpix2

    return x, y


def xy2skypv(x, y, crpix1, crpix2, crval1, crval2, cd, pv, nord):
    """
    Transforms from pixel coordinates to celestial coordinates taking
    non-linear distortion into account with the World Coordinate System
    FITS keywords as used in MegaPipe.

    For the inverse operation see sky2xy

    Reference material:
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/megapipe/docs/CD_PV_keywords.pdf

    Args:
      x, y: float
        Input pixel coordinate
      crpix1: float
        Tangent point x, pixels
      crpix2: float
        Tangent point y, pixels
      crval1: float
        Tangent point RA, degrees
      crval2:
        Tangent point Dec, degrees
      cd: 2d array
        Expresses the scale, the rotation and any possible skew of the image
        with respect to the sky.
      pv: 2d array
      nord: int
        order of the fit

    Returns:
      ra: float
        Right ascension
      dec: float
        Declination
    """
    xp = x - crpix1
    yp = y - crpix2

    # IMPORTANT NOTE: 0-based indexing in Python means indexing for values
    # in cd and pv will be shifted from in the paper.
    x_deg = cd[0][0] * xp + cd[0][1] * yp
    y_deg = cd[1][0] * xp + cd[1][1] * yp

    if nord < 0:
        xi = x
        eta = y
    else:
        xi = pv[0][0]
        eta = pv[1][0]

    if nord >= 1:
        r = numpy.sqrt(x_deg ** 2 + y_deg ** 2)
        xi += pv[0][1] * x_deg + pv[0][2] * y_deg + pv[0][3] * r
        eta += pv[1][1] * y_deg + pv[1][2] * x_deg + pv[1][3] * r

        if nord >= 2:
            x2 = x_deg ** 2
            xy = x_deg * y_deg
            y2 = y_deg ** 2

            xi += pv[0][4] * x2 + pv[0][5] * xy + pv[0][6] * y2
            eta += pv[1][4] * y2 + pv[1][5] * xy + pv[1][6] * x2

            if nord >= 3:
                x3 = x_deg ** 3
                x2y = x2 * y_deg
                xy2 = x_deg * y2
                y3 = y_deg ** 3

                xi += pv[0][7] * x3 + pv[0][8] * x2y + pv[0][9] * xy2 + pv[0][10] * y3
                eta += pv[1][7] * y3 + pv[1][8] * xy2 + pv[1][9] * x2y + pv[1][10] * x3

    xir = xi / PI180
    etar = eta / PI180

    ra0 = crval1 / PI180
    dec0 = crval2 / PI180

    ctan = numpy.tan(dec0)
    ccos = numpy.cos(dec0)
    raoff = numpy.arctan2(xir / ccos, 1 - etar * ctan)
    ra = raoff + ra0
    dec = numpy.arctan(numpy.cos(raoff) / ((1 - (etar * ctan)) / (etar + ctan)))

    ra *= PI180
    ra = numpy.where(ra < 0 , ra + 360, ra)
    ra = numpy.where(ra > 360, ra - 360, ra)
    # f = numpy.int(ra < 0) + 360
    # print f
    # if ra < 0:
    #     ra += 360
    # if ra >= 360:
    #     ra -= 360

    dec *= PI180

    return ra * units.degree, dec * units.degree


def parse_cd(header):
    """
    Parses the CD array from an astropy FITS header.

    Args:
      header: astropy.io.fits.header.Header
        The header containing the CD values.

    Returns:
      cd: 2d array (list(list(float))
        [[CD1_1, CD1_2], [CD2_1, CD2_2]]
    """
    return [[header["CD1_1"], header["CD1_2"]],
            [header["CD2_1"], header["CD2_2"]]]


def parse_pv(header):
    """
    Parses the PV array from an astropy FITS header.

    Args:
      header: astropy.io.fits.header.Header
        The header containing the PV values.

    Returns:
      cd: 2d array (list(list(float))
        [[PV1_0, PV1_1, ... PV1_N], [PV2_0, PV2_1, ... PV2_N]]
        Note that N depends on the order of the fit.  For example, an
        order 3 fit goes up to PV?_10.
    """
    order_fit = parse_order_fit(header)

    def parse_with_base(i):
        key_base = "PV%d_" % i

        pvi_x = [header[key_base + "0"]]

        def parse_range(lower, upper):
            for j in range(lower, upper + 1):
                pvi_x.append(header[key_base + str(j)])

        if order_fit >= 1:
            parse_range(1, 3)

        if order_fit >= 2:
            parse_range(4, 6)

        if order_fit >= 3:
            parse_range(7, 10)

        return pvi_x

    return [parse_with_base(1), parse_with_base(2)]


def parse_order_fit(header):
    """
    Parses the order of the fit for PV terms.

    Args:
      header: astropy.io.fits.header.Header
        The header containing the PV values.

    Returns:
      order_fit: int
        The order of the fit.
    """
    return int(header["NORDFIT"])
