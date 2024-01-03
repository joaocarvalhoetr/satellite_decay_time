import math
from constants import *

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

def density_calculate(altitude):
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

# Example usage
altitude = 0
temperature = temperature_calculate(altitude)
scale_factor = scale_factor_calculate(altitude)
density = density_calculate(altitude)

print(f"Altitude: {altitude} km")
print(f"Temperature: {temperature} K")
print(f"Scale factor: {scale_factor} m")
print(f"Density: {density} kg/m^3")
