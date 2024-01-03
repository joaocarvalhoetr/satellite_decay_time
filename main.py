import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
from atmosphericlayers import *
from statevectorcalculation import *
from extrema import extrema
from scipy.signal import argrelextrema
from scipy.integrate import odeint
from density import atmosphere
from findpeaks import *
from trid_plot import *
from teste_b import *
from ciraModel import *

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
time_vec = []
radious = []

#orbit_tridimensional(a,T,e,RA ,i,w )


# Function to calculate the acceleration
def acceleration(t, y):
    r = y[:3]
    v = y[3:]
    drdt = v

    # Constants
    mu = 398600  # Gravitational parameter in km^3/s^2
    R_earth = 6378  # Earth radius in km
    wE = np.array([0, 0, 7.2921159e-5])  # Earth's angular velocity in rad/s
    CD = 2.2  # Drag coefficient
    A = np.pi / 4 * (1 ** 2)  # Frontal area in m^2
    m = 100  # Mass in kg

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

    # append time and altitude to the lists
    time_vec.append(t)
    radious.append(altitude)

    if altitude <= 0:
        # if the altitude is less than or equal to 0, the integration stops
        print("Termination event triggered.")
        return 0
    return 1

# Integration Settings

t0 = 1e-100
tf = 120 * days

R0, V0 = sv_from_coe(initial_vec, mu)

# initial state vector
y0 = [R0[0][0], R0[0][1], R0[0][2], V0[0][0], V0[0][1], V0[0][2]]


# Number of points to output
nout = 40000

# Integration Time Interval from t0 to tf with nout points
tspan = np.linspace(t0, tf, nout)

# Set error tolerances, initial step size, and termination event:
options = {'rtol': 1e-15, 'atol': 1e-15}

event.terminal = True

# Print initial altitude
initial_altitude = np.linalg.norm(y0[:3]) - R_earth
print("Initial Altitude:", initial_altitude)

# Call the ODE solver

#sol= solve_ivp(acceleration, [t0, tf], y0, method='RK45', t_eval=tspan, events=event)
sol = solve_ivp(acceleration, (t0, tf), y0, t_eval=tspan, events=event, method='DOP853')

# pick all the y values from sol, compute for each y the norm, and subtract the earth radius
altitude = np.linalg.norm(sol.y[:3], axis=0) - R_earth

time = sol.t

[maxima, minima] = find_local_extrema(altitude, time)

# Example usage
maxima_x = time[maxima]
maxima_y = altitude[maxima]

minima_x = time[minima]
minima_y = altitude[minima]

# Fit curves to maxima and minima
degree = 15  # You can adjust the degree of the polynomial
maxima_coefficients, maxima_fitted_values = fit_curve(maxima_x, maxima_y, degree)
minima_coefficients, minima_fitted_values = fit_curve(minima_x, minima_y, degree)

# Plotting maxima with fitted values
plt.plot(time[maxima]/8640, maxima_fitted_values, label='Maxima', linestyle='-', color='red')
# Plotting minima with fitted values
plt.plot(time[minima]/8640, minima_fitted_values, label='Minima', linestyle='-', color='green')

# Adding labels and legend
plt.xlabel('Time (days)')  # Replace 'Time' with the actual label for the x-axis
plt.ylabel('Fitted Values')  # Replace 'Fitted Values' with the actual label for the y-axis
plt.legend()

# Display the plot
plt.show()

