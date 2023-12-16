import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import os

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
r = np.array([0, 0, 0])  # Initialize position vector
v = np.array([0, 0, 0])  # Initialize velocity vector

# Function to calculate spacecraft acceleration
def rates(t, f):
    R = f[:3]  # Position vector (km)
    r = np.linalg.norm(R)  # Distance from earth’s center (km)
    alt = r - R_earth  # Altitude (km)
    V = f[3:]  # Velocity vector (km/s)

    V = np.squeeze(V)
    R = np.squeeze(R)

    Vrel = V - np.cross(wE, R)  # Velocity relative to the atmosphere (km/s)

    vrel = np.linalg.norm(Vrel)  # Speed relative to the atmosphere (km/s)
    uv = Vrel / vrel  # Relative velocity unit vector
    ap = -CD * A / m * 1.225 * (1000 * vrel) ** 2 / 2 * uv  # Acceleration due to drag (m/s^2)
    a0 = -mu * R / r ** 3  # Gravitational acceleration (km/s^2)
    a = a0 + ap / 1000  # Total acceleration (km/s^2)
    return np.concatenate([V, a])

# Initial conditions
y0 = np.concatenate([r, v]).ravel()  # Flatten the array

# Integration time interval
t0, tf = 0, 120 * 86400  # Initial and final times in seconds
nout = 40000  # Number of solution points to output
t_span = np.linspace(t0, tf, nout)

# Set error tolerances, initial step size, and termination event
options = {
    'rtol': 1e-8,
    'atol': 1e-8,
    'first_step': (tf - t0) / 10000,
}

# Initialize lists to store apogee and perigee altitudes and times
apogee_altitudes = []
perigee_altitudes = []
apogee_times = []
perigee_times = []

def altitude_event(t, y):
    R = y[:3]
    
    # Altitude
    altitude = np.linalg.norm(R) - R_earth
    #altitude event print
    print(altitude)
    # Check for apogee and perigee
    if altitude >= max(apogee_altitudes, default=float('-inf')):
        apogee_altitudes.append(altitude)
        apogee_times.append(t)

    if altitude <= min(perigee_altitudes, default=float('inf')):
        perigee_altitudes.append(altitude)
        perigee_times.append(t)

    return altitude

altitude_event.terminal = True

# Integration using solve_ivp
sol = solve_ivp(rates, [t0, tf], y0, t_eval=t_span, method='RK45', vectorized=True, events=altitude_event, **options)


# Plot altitudes within the specified time range
# Extract the locally extreme altitudes within the specified time range
time_range_indices = (np.array(sol.t) >= t0) & (np.array(sol.t) <= tf)
if any(time_range_indices):
    altitude_range = np.linalg.norm(np.array(sol.y)[:3, time_range_indices], axis=0) - RE
    time_range = np.array(sol.t)[time_range_indices]

    # Plot altitudes within the specified time range
    plt.plot(time_range / 86400, altitude_range, 'k', linewidth=2, label='Altitude')
    plt.scatter(np.array(apogee_times) / 86400, np.array(apogee_altitudes), color='r', marker='o', label='Apogee')
    plt.scatter(np.array(perigee_times) / 86400, np.array(perigee_altitudes), color='b', marker='o', label='Perigee')
    plt.grid(True)
    plt.xlabel('Time (days)')
    plt.ylabel('Altitude (km)')
    plt.legend()
    plt.show()
else:
    print("No altitude events within the specified time range.")

# Remove the file atmosphericlayers_output.csv if it exists
if os.path.exists('atmosphericlayers_output.csv'):
    os.remove('atmosphericlayers_output.csv')
else:
    print("The file atmosphericlayers_output.csv does not exist.")
