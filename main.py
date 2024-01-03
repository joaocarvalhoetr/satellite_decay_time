import math
from constants import *
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
from statevectorcalculation import *
from extrema import extrema
from scipy.signal import argrelextrema
from scipy.integrate import odeint
from density import atmosphere
from findpeaks import *
from trid_plot import *
from teste_b import *
from ciraModel import *

# Read from data.txt, CD, m, A, R_p, R_a, RA, i, w, TA
# The structure of the lines in the txt is for example for Cd "Cd: 2.2"
# After reading the values, compute  e, a, h, T

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
time_vec = []
radious = []

#orbit_tridimensional(a,T,e,RA ,i,w )

# Function to calculate the velocity and acceleration
def vel_and_acceleration(t, y):
    r = y[:3]
    v = y[3:]
    drdt = v

    # Calculate atmospheric velocity
    v_atmosphere = np.cross(wE, r)

    # Calculate relative velocity
    v_rel = v - v_atmosphere

    # Unitary vector
    vrel_unit = v_rel / np.linalg.norm(v_rel)

    # Compute the absolute value of the relative velocity
    vrel_abs = np.linalg.norm(v_rel)

    # Calculate density
    density = atmosphere2(np.linalg.norm(r) - R_earth)  # Adjust this value or use a function to compute density based on altitude

    #print(density)

    # Calculate drag acceleration
    P = -CD *A/m * density * (1000*vrel_abs)**2/2 * vrel_unit

    # Calculate gravitational acceleration
    a0 = (-mu / float(np.linalg.norm(r))**3) * np.array(r)

    # Calculate the total acceleration
    dvdt = a0 + P/1000   

    return np.concatenate((drdt, dvdt))

# Function to detect the end of the integration
def event(t, y):
    altitude = np.linalg.norm(y[:3]) - R_earth
    
    #altitude triggered


    if altitude <= 0:
        # if the altitude is less than or equal to 0, the integration stops
        print("Termination event triggered.")
        return 0
    return 1

# Integration Settings

t0 = 0
tf = 150 * days

R0, V0 = sv_from_coe(initial_vec, mu)

# initial state vector
y0 = [R0[0][0], R0[0][1], R0[0][2], V0[0][0], V0[0][1], V0[0][2]]


# Number of points to output
nout = 40000

# Integration Time Interval from t0 to tf with nout points
tspan = np.linspace(t0, tf, nout)

# Set error tolerances, initial step size, and termination event:
options = {'rtol': 1e-7, 'atol': 1e-7}

event.terminal = True

# Print initial altitude
initial_altitude = np.linalg.norm(y0[:3]) - R_earth
print("Initial Altitude:", initial_altitude)

# Call the ODE solver
# For each iteration, call the events function to check if the integration should stop or not.

sol = solve_ivp(vel_and_acceleration, (t0, tf), y0, events=event, method='DOP853')

# Pick all the y values from sol, compute for each y the norm, and subtract the earth radius.

altitude = np.linalg.norm(sol.y[:3], axis=0) - R_earth

# Time for each iteration
time = sol.t *10

# Compute the extrema for the data aquired

[maxima, minima] = find_local_extrema(altitude, time)

# Organize the data for the maxima and minima

maxima_x = time[maxima]
maxima_y = altitude[maxima]

minima_x = time[minima]
minima_y = altitude[minima]

# Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs

degree = 15  # Degree of the polynomial curve to fit the data
maxima_coefficients, maxima_fitted_values = fit_curve(maxima_x, maxima_y, degree)
minima_coefficients, minima_fitted_values = fit_curve(minima_x, minima_y, degree)

# Plotting maxima with fitted values
plt.plot(time[maxima]/86400, maxima_fitted_values, label='Maxima', linestyle='-', color='red')

# Plotting minima with fitted values
plt.plot(time[minima]/86400, minima_fitted_values, label='Minima', linestyle='-', color='green')

# Adding labels and legend
plt.xlabel('Time (days)')  # Replace 'Time' with the actual label for the x-axis
plt.ylabel('Fitted Values')  # Replace 'Fitted Values' with the actual label for the y-axis
plt.legend()

# Display the plot
plt.show()

