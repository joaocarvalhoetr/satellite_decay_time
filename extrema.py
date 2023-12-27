import numpy as np
from scipy.signal import argrelextrema

def extrema(x):
    # Find local minima and maxima indices
    maxima_idx = argrelextrema(x, np.greater)[0]
    minima_idx = argrelextrema(x, np.less)[0]

    # Get corresponding values
    xmax = x[maxima_idx]
    imax = maxima_idx
    xmin = x[minima_idx]
    imin = minima_idx

    # Sort in descending order
    descending_order = np.argsort(xmax)[::-1]
    xmax = xmax[descending_order]
    imax = imax[descending_order]

    ascending_order = np.argsort(xmin)
    xmin = xmin[ascending_order]
    imin = imin[ascending_order]

    return xmax, imax, xmin, imin