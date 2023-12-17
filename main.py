import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import os
from atmosphericlayers import *
from statevectorcalculation import *

# Constants
RE = 6371.0  # Earth's radius in km
mu = 398600.4418  # Earth's gravitational parameter in km^3/s^2
wE = [0, 0, 7.2921159e-5]  # Earth's angular velocity in rad/s
deg = math.pi / 180.0  # Conversion factor from degrees to radians
R_earth = 6371.0  # Earth's radius in km

# Read from data.txt, Cd(1st line), e(2nd line), and R_a(3rd line)
# The structure of the lines in the txt is for example for Cd "Cd: 2.2"
with open('data.txt') as f:
    lines = f.readlines()

    # Drag Coefficient
    CD = float(lines[0].split()[1])

    # Mass of the satellite
    m = float(lines[1].split()[1])

    # Area of the satellite
    A = float(lines[2].split()[1])

    # Eccentricity
    R_p = float(lines[3].split()[1]) + R_earth

    # Apoapsis radius
    R_a = float(lines[4].split()[1]) + R_earth

    # Right ascension of the ascending node
    RA = float(lines[5].split()[1]) * deg

    # Inclination
    i = float(lines[6].split()[1]) * deg

    # Argument of perigee
    w = float(lines[7].split()[1]) * deg

    # True Anomaly
    TA = float(lines[8].split()[1]) * deg

    # Eccentric anomaly
    e = (R_a - R_p) / (R_a + R_p)

    # Semi-major axis
    a = (R_a + R_p) / 2

    # Angular momentum
    h = math.sqrt(mu * a * (1 - e ** 2))

    # Period
    T = 2 * math.pi * math.sqrt(a ** 3 / mu)

# Calculate initial state vector
initial_vec = [h, e, RA, i, w, TA]

# Function to calculate the acceleration
def acceleration(t, y):
    r = y[:3]
    v = y[3:]
    
    # Calculate atmospherical velocity
    wE = [0, 0, 7.2921159e-5]  # Earth's angular velocity in rad/s
    v_atmosphere = np.cross(wE, r)

    # Calculate relative velocity
    v_rel = v - v_atmosphere

    # unitary vector
    vrel_unit = v_rel / np.linalg.norm(v_rel)

    # compute the absolute value of the relative velocity
    vrel_abs = np.linalg.norm(v_rel)

    # Calculate P
    P = -0.5 * density_calculate(np.linalg.norm(r)-R_earth) * CD * A / m * np.linalg.norm(v_rel) * vrel_unit

    # Calculate gravitational acceleration
    a0 = -mu / np.linalg.norm(r) ** 3 * r

    # Calculate the acceleration
    a = a0 + P + (-mu / np.linalg.norm(r) ** 3) * r

    return np.concatenate((v, a))

def event(t, y):
    # Print altitude during integration
    altitude = np.linalg.norm(y[:3]) - R_earth
    print("Altitude:", altitude)

    # return the difference between altitude and 0
    return altitude

# Integration Settings
t0 = 0
tf = 120 * days

R0, V0 = sv_from_coe(initial_vec, mu)

print("Initial state vector: ", R0, V0)

# initial state vector
y0 = [R0[0][0], R0[0][1], R0[0][2], V0[0][0], V0[0][1], V0[0][2]]

# Number of points to output
nout = 40000

# Integration Time Interval from t0 to tf with nout points
tspan = np.linspace(t0, tf, nout)

# Set error tolerances, initial step size, and termination event:
options = {'rtol': 1e-8, 'atol': 1e-8, 'h0': 1/10000}

# Call the ODE solver
sol = solve_ivp(acceleration, [t0, tf], y0, method='RK45', t_eval=tspan, events=event, **options)

# Extract locally extreme altitudes
altitude = np.sqrt(np.sum(sol.y[:3, :] ** 2, axis=0)) - R_earth  # Altitude at each time

# Find local extrema
maxima_indices = np.r_[True, altitude[1:] > altitude[:-1]] & np.r_[altitude[:-1] > altitude[1:], True]
minima_indices = np.r_[True, altitude[1:] < altitude[:-1]] & np.r_[altitude[:-1] < altitude[1:], True]

# Get times and altitudes of maxima and minima
maxima_times = sol.t[maxima_indices]
maxima_altitudes = altitude[maxima_indices]
minima_times = sol.t[minima_indices]
minima_altitudes = altitude[minima_indices]

# Create arrays for maxima and minima
maxima = np.column_stack((maxima_times, maxima_altitudes))
minima = np.column_stack((minima_times, minima_altitudes))

# Sort maxima and minima by time
apogee = maxima[maxima[:, 0].argsort()]
perigee = minima[minima[:, 0].argsort()]

# Plot perigee and apogee history on the same figure
plt.figure(1)

# Set NaN for the first value in apogee (assuming it's a 2D array)
apogee[0, 1] = np.nan

# Plot apogee in blue
plt.plot(apogee[:, 0] / days, apogee[:, 1], 'b', linewidth=2, label='Apogee')

# Plot perigee in red
plt.plot(perigee[:, 0] / days, perigee[:, 1], 'r', linewidth=2, label='Perigee')

# Set up the plot with grid, labels, and limits
plt.grid(True)
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.xlabel('Time (days)')
plt.ylabel('Altitude (km)')
plt.ylim([0, 1000])

# Show legend
plt.legend()

# Display the plot
plt.show()


# Remove the file atmosphericlayers_output.csv if it exists
if os.path.exists('atmosphericlayers_output.csv'):
    os.remove('atmosphericlayers_output.csv')
else:
    print("The file atmosphericlayers_output.csv does not exist.")
