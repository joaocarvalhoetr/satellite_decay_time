import numpy as np

def sv_from_coe(coe, mu):
    h, e, RA, incl, w, TA = coe

    rp = (h**2 / mu) * (1 / (1 + e * np.cos(TA))) * np.array([np.cos(TA), np.sin(TA), 0])
    vp = (mu / h) * np.array([-np.sin(TA), e + np.cos(TA), 0])

    R3_W = np.array([[np.cos(RA), np.sin(RA), 0],
                     [-np.sin(RA), np.cos(RA), 0],
                     [0, 0, 1]])

    R1_i = np.array([[1, 0, 0],
                     [0, np.cos(incl), np.sin(incl)],
                     [0, -np.sin(incl), np.cos(incl)]])

    R3_w = np.array([[np.cos(w), np.sin(w), 0],
                     [-np.sin(w), np.cos(w), 0],
                     [0, 0, 1]])

    Q_pX = np.transpose(np.dot(R3_w, np.dot(R1_i, R3_W)))

    r = np.dot(Q_pX, rp)
    v = np.dot(Q_pX, vp)

    r = r.reshape((1, 3))
    v = v.reshape((1, 3))

    return r, v

# Example data
deg = np.pi / 180
mu = 398600
h = 80000
e = 1.4
RA = 40 * deg
incl = 30 * deg
w = 60 * deg
TA = 30 * deg

coe = [h, e, RA, incl, w, TA]


