import numpy as np


def pos_cart(a: np.array, e: np.array, inc: np.array, Omega: np.array, omega: np.array, M: np.array) -> np.array:
    """
    This routine transforms delaunay variables into cartesian variables, positions only.

    Parameters
    a: semi-major axis (AU)
    e: eccentricity
    inc: inclination (radians)
    Omega: longitude of the ascending node (radians)
    omega: argument of the pericentre (radians)
    M: mean anomaly (radians)

    """
    M = M % (2 * np.pi)
    signe = np.sign(a)
    cos_i = np.cos(inc)
    sin_i = np.sqrt(1 - cos_i ** 2)
    delau1 = M
    delau2 = np.cos(omega)
    delau3 = np.sin(omega)
    delau4 = np.cos(Omega)
    delau5 = np.sin(Omega)
    delau6 = signe * np.sqrt(a * signe)
    delau7 = np.abs(delau6) * np.sqrt((1 - e ** 2) * signe)

    E = compute_E(e, M)

    cos_E = np.cos(E)
    sin_E = np.sin(E)
    q_vec = np.zeros_like(delau3, shape=(2, delau6.shape[0]))
    q_vec[0,:] = delau6 ** 2 * (cos_E - e)
    q_vec[1,:] = delau7 * delau6 * sin_E

    mat = np.zeros_like(delau3, shape=(3, 2, delau1.shape[0]))
    d53 = delau5 * delau3
    d42 = delau4 * delau2
    d52 = delau5 * delau2
    d43 = delau4 * delau3
    # Rotation matrix must match F95/elemutils.f95 pos_cart (signs matter).
    mat[0, 0, :] = d42 - cos_i * d53
    mat[0, 1, :] = -d43 - cos_i * d52
    mat[1, 0, :] = d52 + cos_i * d43
    mat[1, 1, :] = -d53 + cos_i * d42
    mat[2, 0, :] = sin_i * delau3
    mat[2, 1, :] = sin_i * delau2

    # Cartesian coordinates
    p = (mat * q_vec).sum(axis=1)
    # p = np.zeros(3)
    # p[0] = mat[0, 0] * q_vec[0] + mat[0, 1] * q_vec[1]
    # p[1] = mat[1, 0] * q_vec[0] + mat[1, 1] * q_vec[1]
    # p[2] = mat[2, 0] * q_vec[0] + mat[2, 1] * q_vec[1]

    return p


def compute_E(e: np.array, M: np.array) -> np.array:
    """
    Compute the eccentric anomaly, E, from the mean anomaly, M, and the eccentricity, e.
    """
    E = M + 0.85 * np.sign(np.sin(M)) * e
    i = 0
    while True:
        sin_e = e * np.sin(E)
        f = E - sin_e - M
        iterate = np.fabs(f) > 1e-14
        if iterate.sum() > 0:
            cos_e = e[iterate] * np.cos(E[iterate])
            fp = 1 - cos_e
            fpp = sin_e[iterate]
            fppp = cos_e
            de = -f[iterate] / fp
            de = -f[iterate] / (fp + de * fpp / 2)
            de = -f[iterate] / (fp + de * fpp / 2 + de * de * fppp / 6)
            E[iterate] = E[iterate] + de
            i += 1
            if i < 20:
                continue
            raise ValueError('POS_CART: No convergence after {i} iterations')
        break
    return E


if __name__ == '__main__':
    a = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    e = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    inc = np.array([0.0, np.pi/4.0, np.pi/3, np.pi/2, np.pi])
    omega = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    Omega = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    M = np.array([0.0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])

    print(pos_cart(a, e, inc, Omega, omega, M))