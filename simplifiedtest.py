import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from atmosphericlayers import *

Cd = 2.2
m = 100
A = 1
altitude = 650

R0 = altitude + R_earth

def f(t, y):
    x, y, vx, vy = y
    r = math.sqrt(x**2 + y**2)
    vrel = math.sqrt(vx**2 + vy**2) - omega_earth*r

    ax = -mu*x/r**3 - 0.5* density_calculate(r- R_earth)* vrel**2 * Cd * A/m
    ay = -mu*y/r**3 - 0.5* density_calculate(r- R_earth)* vy**2 * Cd * A/m

    return [vx, vy, ax, ay]

def main():
    x0 = 0
    y0 = R0
    vx0 = math.sqrt(mu/y0)
    vy0 = 0

    t0 = 0
    tf = 100000
    dt = 0.1
    t = np.arange(t0, tf, dt)

    sol = solve_ivp(f, [t0, tf], [x0, y0, vx0, vy0], t_eval=t)

    r = np.sqrt(sol.y[0]**2 + sol.y[1]**2) - R_earth

    # Identify indices where altitude is negative
    neg_altitude_indices = np.where(r < 0)[0]

    # Exclude negative altitude points from the plot
    plt.plot(np.delete(sol.t, neg_altitude_indices), np.delete(r, neg_altitude_indices))

    # Grid on the plot
    plt.grid()

    # Label the axes
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (km)')

    plt.show()

main()
