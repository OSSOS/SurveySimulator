import numpy as np
import pylab as plt


def el2rv(mu, a, e, inc, capom, om, f):
    """
    Takes orbital elements and calculates position
    and velocity vectors.
    Input should be in ***RADIANS***
    Based on C code written by Brett Gladman,
    translated to python by Sam Lawler, summer 2014
    """

    prec = 1.0e-13  # user can change this if more precision needed (just runs slower)

    # compute the unit vector
    u = om + f
    xhat = np.cos(u) * np.cos(capom) - np.cos(inc) * np.sin(capom) * np.sin(u)
    yhat = np.cos(u) * np.sin(capom) + np.cos(inc) * np.cos(capom) * np.sin(u)
    zhat = np.sin(inc) * np.sin(u)

    # compute the angular momentum vector (unit vector)
    hx = np.sin(capom) * np.sin(inc)
    hy = -np.cos(capom) * np.sin(inc)
    hz = np.cos(inc)

    # assuming not parabolic, here the magnitudes of the vectors
    r = a * (1.0 - e * e) / (1.0 + e * np.cos(f))
    h = (mu * a * (1.0 - e * e)) ** 0.5

    # position vectors
    x = r * xhat
    y = r * yhat
    z = r * zhat

    # compute components of vector theta hat
    thx = hy * zhat - hz * yhat
    thy = hz * xhat - hx * zhat
    thz = hx * yhat - hy * xhat

    # obtain the velocity vector's components and calculate v
    thdot = h / (r * r)
    rdot = e * mu * np.sin(f) / h

    vx = r * thdot * thx + rdot * xhat
    vy = r * thdot * thy + rdot * yhat
    vz = r * thdot * thz + rdot * zhat

    return x, y, z


def mtof(e, m):
    """converts mean anomaly to true anomaly for arrays
  Input should be in ***RADIANS*** """
    # first calculate eccentric anomaly (bigE)
    f = np.zeros(len(e))
    for ind in np.arange(0, len(e)):
        n = 0.
        delta = 1000.
        big_e = m[ind] - e[ind] * np.sin(m[ind])
        while n < 1.e4 and delta > 1.e-6:
            f1 = big_e - e[ind] * np.sin(big_e) - m[ind]
            fp = 1.0 - e[ind] * np.cos(big_e)
            delta = -f1 / fp
            big_e = big_e + delta
            n = n + 1
        f[ind] = 2. * np.arctan(((1. + e[ind]) / (1. - e[ind])) ** 0.5 * np.tan(big_e / 2.))
    return f


#####################
da, de, dinc, dnode, dargperi, dmanom = np.genfromtxt("drawn.dat", usecols=(0, 1, 2, 3, 4, 5), unpack=True)
ta, te, tinc, tnode, targperi, tmanom = np.genfromtxt("tracked.dat", usecols=(0, 1, 2, 3, 4, 5), unpack=True)
df = mtof(de, dmanom / 180. * np.pi)
tf = mtof(te, tmanom / 180. * np.pi)

dx = []
dy = []
tx = []
ty = []
for i in np.arange(0, len(da)):
    dx1, dy1, dz1 = el2rv(1.0, da[i], de[i], dinc[i] / 180. * np.pi, dnode[i] / 180. * np.pi,
                          dargperi[i] / 180. * np.pi, df[i])
    # print (x1**2.+y1**2.)**0.5
    dx = np.append(dx, dx1)
    dy = np.append(dy, dy1)
for i in np.arange(0, len(ta)):
    tx1, ty1, tz1 = el2rv(1.0, ta[i], te[i], tinc[i] / 180. * np.pi, tnode[i] / 180. * np.pi,
                          targperi[i] / 180. * np.pi, tf[i])
    # print (x1**2.+y1**2.)**0.5
    tx = np.append(tx, tx1)
    ty = np.append(ty, ty1)

# Neptune's postion
lambdaN = 333.5078 * np.pi / 180.
xN = 30.1 * np.cos(lambdaN)
yN = 30.1 * np.sin(lambdaN)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
plt.plot(dx, dy, '.k', markersize=1.)
plt.plot(tx, ty, '.r', markersize=3.)
plt.plot(xN, yN, 'ob', mec='k')
plt.plot(0, 0, 'oy', mec='k')
circle1 = plt.Circle((0, 0), 30., fill=False, ec='b', ls=':')
circle2 = plt.Circle((0, 0), 40., fill=False, ec='b', ls=':')
circle3 = plt.Circle((0, 0), 50., fill=False, ec='b', ls=':')
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
plt.text(0, 31, "30 AU", horizontalalignment='center', color='b')
plt.text(0, 41, "40 AU", horizontalalignment='center', color='b')
plt.text(0, 51, "50 AU", horizontalalignment='center', color='b')
plt.axis('off')

plt.savefig('example.png')
plt.show()
