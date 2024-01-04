import numpy as np

def find_local_extrema(altitude, time):
    # Calculate the derivative
    derivative = np.diff(altitude)

    # Identify zero-crossings
    zero_crossings = np.where(np.diff(np.sign(derivative)))[0]

    # Refine maxima and minima based on the second derivative
    second_derivative = np.diff(derivative)
    maxima = zero_crossings[second_derivative[zero_crossings - 1] < 0]
    minima = zero_crossings[second_derivative[zero_crossings - 1] > 0]

    # Handle edge cases
    if altitude[-1] > altitude[-2]:
        maxima = np.append(maxima, len(altitude) - 1)
    elif altitude[-1] < altitude[-2]:
        minima = np.append(minima, len(altitude) - 1)

    return maxima, minima


def fit_curve(x, y, degree):
    # Fit a polynomial curve using least squares
    coefficients = np.polyfit(x, y, degree)
    
    # Generate fitted values
    fitted_values = np.polyval(coefficients, x)

    return coefficients, fitted_values
