import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import the 3D plotting toolkit
from matplotlib.animation import FuncAnimation
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
time_vec = []
radious = []


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

# Function to detect the end of the integration
def event(t, y):
    

    # Print altitude during integration
    altitude = np.linalg.norm(y[:3]) - R_earth

    #append time and altitude to the lists
    time_vec.append(t)
    radious.append(altitude)

    if altitude <= 0:
        # if the altitude is less than or equal to 0, the integration stops
        return 0
    return 1


# Integration Settings
t0 = 0
tf = 500000 * days

R0, V0 = sv_from_coe(initial_vec, mu)

print("Initial state vector: ", R0, V0)

# initial state vector
y0 = [R0[0][0], R0[0][1], R0[0][2], V0[0][0], V0[0][1], V0[0][2]]

# Number of points to output
nout = 40000

# Integration Time Interval from t0 to tf with nout points
tspan = np.linspace(t0, tf, nout)

# Set error tolerances, initial step size, and termination event:
options = {'rtol': 1e-10, 'atol': 1e-10}

event.terminal = True

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Call the ODE solver
sol = solve_ivp(acceleration, [t0, tf], y0, method='RK45', t_eval=tspan, events=event, **options, dense_output=True)

def update_plot(i, time_vec, radious):
    ax.clear()
    ax.set_xlim(0, 250)  # Adjust these limits based on your simulation
    ax.set_ylim(0, 250)
    ax.set_zlim(0, 250)
    ax.plot(time_vec[:i], radious[:i], zs=0, zdir='z', label='Satellite Trajectory')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Altitude (km)')
    ax.set_zlabel('Z')

# Animation function
def update(i):
    update_plot(i, time_vec, radious)

# Create animation
animation = FuncAnimation(fig, update_plot, frames=len(sol.t), fargs=(time_vec, radious), interval=50, repeat=False)

plt.show()

#plot radious vs time
plt.plot(time_vec, radious)
plt.xlabel('Time (s)')
plt.ylabel('Altitude (km)')
plt.title('Altitude vs Time')
plt.show()

# Remove the file atmosphericlayers_output.csv if it exists
if os.path.exists('atmosphericlayers_output.csv'):
    os.remove('atmosphericlayers_output.csv')
else:
    print("The file atmosphericlayers_output.csv does not exist.")
