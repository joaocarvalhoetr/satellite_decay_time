import numpy as np
import math
from constants import *

#Student Made

# Given constants
R_air = 287.05  # gas constant for air
g_0 = 9.80665  # sea level gravitational acceleration
r_0 = 1.225 # sea level density

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
    density = r_0 * math.exp(-(altitude*1000)/ scale_factor)
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
     
#Gram

#Geometric Altitudes
h = np.array([145, 150, 155, 160, 170, 180, 190, 200, 210, 220,
 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 
370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620,
 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900])

 # Corresponding densities (kg/m^3) for T_Exosphere = 2400k
r = np.array([3.39e-9 , 2.63e-9, 2.08e-9, 1.69e-9, 1.17e-9, 8.48e-10, 6.44e-10, 5.05e-10,
4.06e-10, 3.34e-10, 2.79e-10, 2.36e-10, 2.02e-10, 1.75e-10, 1.52e-10, 1.33e-10, 1.17e-10, 1.04e-10,
9.19e-11, 8.19e-11, 7.32e-11, 6.55e-11, 5.89e-11, 5.30e-11, 4.78e-11, 4.32e-11, 3.91e-11, 3.54e-11,
2.92e-11, 2.43e-11, 2.02e-11, 1.70e-11, 1.43e-11, 1.21e-11, 1.02e-11, 8.70e-12, 7.43e-12,
7.43e-12, 6.36e-12, 5.46e-12, 4.70e-12, 4.06e-12, 3.51e-12, 3.04e-12, 2.64e-12, 2.30e-12, 2.01e-12, 1.75e-12,1.54e-12,
1.35e-12, 1-18e-12, 1-04e-12, 9.16e-13, 8.08e-13])


