import math
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
from atmosphericlayers import *

Cd = 2.2
m = 100
A =1
altitude = 1500

R0 = altitude + R_earth

def f(t, y):

    # y = [x, y, vx, vy]
    # y' = [vx, vy, ax, ay]

    x, y, vx, vy = y
    r = math.sqrt(x**2 + y**2)
    vrel = math.sqrt(vx**2 + vy**2) - omega_earth*r

    ax = -mu*x/r**3 - 0.5* density_calculate(r- R_earth)* vrel**2 * Cd * A/m
    ay = -mu*y/r**3 - 0.5* density_calculate(r- R_earth)* vy**2 * Cd * A/m

    return [vx, vy, ax, ay]

def main():

    # Initial conditions
    x0 = 0
    y0 = R0
    vx0 = math.sqrt(mu/y0)
    vy0 = 0

    # Time
    t0 = 0
    tf = 100000
    dt = 0.1
    t = np.arange(t0, tf, dt)

    # Solve
    sol = solve_ivp(f, [t0, tf], [x0, y0, vx0, vy0], t_eval=t)

    print(sol.t)
    print(sol.y)

    r = np.sqrt(sol.y[0]**2 + sol.y[1]**2) - R_earth

    # Delete all t values and y values for negative altitudes
    for i in range(len(r)):
        if r[i] < 0:
            r = np.delete(r, i)
            sol.t = np.delete(sol.t, i)


    # Plot
    plt.plot(sol.t, r)
    plt.show()


main()
