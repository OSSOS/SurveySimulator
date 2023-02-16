from matplotlib import pyplot
import numpy
import math

def merr(mag, params):
    """
    merr as a function of magnitude given some parameters

    mag:  the magnitude at which the uncertainty estimate is needed
    params:  array of 6 values describing the magnitude uncertainty funciton.
    """

    sqrt = numpy.sqrt
    mag_th = mag
    magerr_bright = params[0]
    magerr_mid =    params[0]*10.0**(params[1]*(mag_th   -21.0))
    magerr_faint =  params[0]*10.0**(params[1]*(params[2]-21.0)) - (mag_th - params[2])*params[3]
    magerr_faint[magerr_faint < 0] = 0

    tmp = numpy.random.random_sample(len(mag_th))
    A = math.sqrt(6.0)*(numpy.sqrt(2*tmp) - 1)
    B = math.sqrt(6.0)*(1 - numpy.sqrt(2*(1-tmp)))
    tmp[tmp <= 0.5] = A[tmp <= 0.5]
    tmp[tmp > 0.5] = B[tmp > 0.5]
    
    magerr = magerr_bright*(mag <= 21.0) + magerr_mid*((mag > 21.0) & ( mag <= params[2])) + magerr_faint*(mag > params[2])
    mag = mag_th + magerr*tmp
    mag[mag_th > params[4]] += (mag_th[mag_th > params[4]] - params[4])*params[5]
    return mag


if __name__ == '__main__':
    params = [0.00,0.00,26.7,-0.00,26.5,-0.0]
    mags_th = numpy.arange(21,27,0.01)
    mags = merr(mags_th, params)

    pyplot.plot(mags_th, mags_th-mags, 'x')
    pyplot.xlabel('mag')
    pyplot.ylabel('mag_err')
    pyplot.show()
