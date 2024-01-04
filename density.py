import numpy as np
import math
from constants import *

#Student Model using Scale Factors

# Temperature data from the initial table
T = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946,
     210.02, 257, 349.49, 892.79, 1022.2, 1103.4, 1205.4, 1322.3,
     1432.1, 1487.4, 1506.1, 1506.1, 1507.6]

# Lapse rates data
L_z = [-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 1.6481, 5.0, 10.0,
       20.0, 15.0, 10.0, 7.0, 5.0, 4.0, 3.3, 2.6, 1.7, 1.1, 0]

# Geometric Altitudes data in km
Z_0 = [0.0, 11.0102, 20.0631, 32.1619, 47.3501, 51.4125, 71.802,
       86.0, 210.65, 260.65, 360.65, 960.65, 1110.65, 1210.65, 1350.65,
       1550.65, 1830.65, 2160.65, 2420.65, 2590.65, 2700.65]

def temperature_calculate(altitude):
    # Locate the layer the given altitude is in
    for i in range(len(Z_0)):
        if Z_0[i] > altitude:
            i -= 1
            break

    i = max(0, i)  # Ensure the index is not negative
    # Compute the temperature for the given altitude
    temperature = T[i] + L_z[i] * (altitude - Z_0[i])
    return temperature

def scale_factor_calculate(altitude):
    # Compute the temperature for the given altitude
    temperature = temperature_calculate(altitude)
    # Compute the scale factor for the given altitude
    scale_factor = (R_air * temperature) / g_0
    return scale_factor

def student_model(altitude):
    # Compute the scale factor for the given altitude
    scale_factor = scale_factor_calculate(altitude)

    # Locate the layer the given altitude is in
    for i in range(len(Z_0)):
        if Z_0[i] > altitude:
            i -= 1
            break

    i = max(0, i)  # Ensure the index is not negative
    # Compute the density for the given altitude
    density = rho_0 * math.exp(-(altitude*1000)/ scale_factor)
    return density


#US1976

def density_US1976(z):
    # Geometric altitudes (km):
    h = np.array([0, 25, 30, 40, 50, 60, 70,
                  80, 90, 100, 110, 120, 130, 140,
                  150, 180, 200, 250, 300, 350, 400,
                  450, 500, 600, 700, 800, 900, 1000])

    # Corresponding densities (kg/m^3) from US Standard Atmosphere 1976:
    r = np.array([1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5,
                  1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9,
                  2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12,
                  1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15])

    # Scale heights (km):
    H = [7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661,
                  5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322,
                  21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019,
                  60.980, 65.654, 76.377, 100.587, 147.203, 208.020]

    # Handle altitudes outside of the range:
    if z > 1000:
        z = 1000
    elif z < 0:
        z = 0

    # Determine the interpolation interval:
    for j in range(27):
        if h[j] <= z < h[j+1]:
            i = j
            break
        elif z == 1000:  # Handle the case when z is exactly 1000
            i = 26

    # Exponential interpolation:
    density = r[i] * np.exp(-(z - h[i]) / H[i])

    return density
    

class MyError(Exception):
    pass


# Jachia+Gram

def density_jacchia(z, T_exo):
    # Geometric Altitudes
    h_ja = np.array([145, 150, 155, 160, 170, 180, 190, 200, 210, 220,
                     230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 
                     370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620,
                     640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 
                     940, 960, 980, 1000])

    h_gram = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 
                       44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 
                       94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 
                       132, 134, 136, 138, 140])

    # Corresponding densities (kg/m^3) for T_Exosphere = 2400k
    r_ja_2400 = np.array([3.39e-9, 2.63e-9, 2.08e-9, 1.69e-9, 1.17e-9, 8.48e-10, 6.44e-10, 5.05e-10,
                          4.06e-10, 3.34e-10, 2.79e-10, 2.36e-10, 2.02e-10, 1.75e-10, 1.52e-10, 1.33e-10, 
                          1.17e-10, 1.04e-10, 9.19e-11, 8.19e-11, 7.32e-11, 6.55e-11, 5.89e-11, 5.30e-11, 
                          4.78e-11, 4.32e-11, 3.91e-11, 3.54e-11, 2.92e-11, 2.43e-11, 2.02e-11, 1.70e-11, 
                          1.43e-11, 1.21e-11, 1.02e-11, 8.70e-12, 7.43e-12, 6.36e-12, 5.46e-12, 4.70e-12, 
                          4.06e-12, 3.51e-12, 3.04e-12, 2.64e-12, 2.30e-12, 2.01e-12, 1.75e-12, 1.54e-12, 
                          1.35e-12, 1-18e-12, 1-04e-12, 9.16e-13, 
                          8.08e-13, 7.13e-13, 6.31e-13, 5.58e-13, 4.95e-13, 4.39e-13])

    # Corresponding densities (kg/m^3) for T_Exosphere = 1500k
    r_ja_1500 = np.array([3.03e-9, 2.31e-9, 1.81e-9, 1.45e-9, 9.81e-10, 6.99e-10, 5.17e-10, 3.94e-10, 
                          3.07e-10, 2.43e-10, 1.95e-10, 1.59e-10, 1.30e-10, 1.08e-10, 8.96e-11, 7.51e-11, 
                          6.32e-11, 5.35e-11, 4.55e-11, 3.89e-11, 3.33e-11, 2.86e-11, 2.47e-11, 2.13e-11,
                          1.85e-11, 1.61e-11, 1.40e-11, 1.23e-11, 9.40e-12, 7.27e-12, 5.66e-12, 4.43e-12, 
                          3.49e-12, 2.76e-12, 2.19e-12, 1.75e-12, 1.40e-12,
                          1.13e-12, 9.07e-13, 7.33e-13, 5.94e-13, 4.83e-13, 3.93e-13, 3.21e-13, 2.63e-13, 
                          2.16e-13, 1.78e-13, 1.47e-13, 1.21e-13, 1.01e-13, 8.39e-14, 7.01e-14, 5.88e-14, 
                          4.94e-14, 4.17e-14, 3.54e-14, 3.01e-14, 2.57e-14])

    # Corresponding densities (kg/m^3) for T_Exosphere = 600k
    r_ja_600 = np.array ([2.07e-9, 1.48e-9, 1.09e-9, 8.15e-10, 4.76e-10, 2.91e-10, 1.83e-10, 1.19e-10, 
                          7.87e-11, 5.32e-11, 3.66e-11, 2.55e-11, 1.80e-11, 1.29e-11, 9.28e-12, 6.74e-12, 
                          4.93e-12, 3.63e-12, 2.68e-12, 1.99e-12, 
                          1.49e-12, 1.11e-12, 8.37e-13, 6.31e-13, 4.77e-13, 3.63e-13, 2.76e-13, 2.12e-13, 
                          1.26e-13, 7.65e-14, 4.79e-14,
                          3.11e-14, 2.10e-14, 1.48e-14, 1.10e-14, 8.46e-15, 6.79e-15, 5.61e-15, 4.75e-15, 
                          4.10e-15, 3.59e-15, 3.17e-15, 2.83e-15, 2.54e-15, 2.28e-15, 2.07e-15, 1.88e-15, 
                          1.71e-15, 1.56e-15, 1.43e-15, 1.31e-15, 1.21e-15, 1.12e-15,
                          1.04e-15, 9.61e-16, 8.94e-16, 8.34e-16, 7.79e-16])

    # Corresponding densities (kg/m^3) for troposphere and stratosphere
    r_gram = np.array([1.225, 9.923e-1, 8.048e-1, 6.525e-1, 5.287e-1, 4.226e-1, 3.291e-1, 
                       2.495e-1, 1.805e-1, 1.335e-1, 9.416e-2, 6.621e-1, 4.698e-2, 3.381e-2, 
                       2.449e-2, 1.783e-2, 1.303e-2,
                       9.581e-3, 7.075e-3, 5.251e-3, 3.924e-3, 2.973e-3, 2.263e-3, 1.742e-3, 
                       1.353e-3, 1.051e-3,
                       8.257e-4, 6.465e-4, 5.059e-4, 3.951e-4, 3.069e-4, 2.365e-4, 1.813e-4, 
                       1.382e-4, 1.049e-4,
                       7.912e-5, 5.905e-5, 4.391e-5, 3.246e-5, 2.389e-5, 1.754e-5, 1.282e-5, 
                       9.341e-6, 6.816e-6, 4.986e-6,
                       3.602e-6, 2.569e-6, 1.814e-6, 1.266e-6, 8.802e-7, 6.114e-7, 4.204e-7, 
                       2.919e-7, 2.033e-7,
                       1.436e-7, 1.035e-7, 7.554e-8, 5.564e-8, 4.309e-8, 3.331e-8, 2.602e-8, 
                       2.057e-8, 1.650e-8,
                       1,342e-8, 1.106e-8, 9.920e-9, 7.773e-9, 6.618e-9, 5.688e-9, 4.929e-9, 
                       4.305e-9 ])

    # Select the correct array for the desired exospheric temperature
    if T_exo == 2400:
        r_ja_selected = r_ja_2400
    elif T_exo == 1500:
        r_ja_selected = r_ja_1500
    elif T_exo == 600:
        r_ja_selected = r_ja_600
    else:
        raise MyError("Exospheric Temperature is not valid.")

    # Combine desired arrays
    r_final = np.concatenate((r_gram, r_ja_selected))
    h_final = np.concatenate((h_gram, h_ja))

    # Handle altitudes outside of the range:
    if z > 1000:
        z = 1000
    elif z < 0:
        z = 0

    # Determine the interpolation interval:
    if z < 180 or z > 500:
        density = np.interp(z, h_final, r_final)
        return density
    else:
        # Select the range of data for fitting
        fit_range_start = 180
        fit_range_end = 500
        mask = (h_final >= fit_range_start) & (h_final <= fit_range_end)
        fit_h = h_final[mask]
        fit_r = r_final[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 10

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)
        
        density = poly_function(z)
        return density

print(density_jacchia(250, 1500))