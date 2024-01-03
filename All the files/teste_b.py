import numpy as np
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
def orbit_tridimensional(a,per,e,Omega,i,w):
    #orbit = pyasl.KeplerEllipse(a=0.62, per=1, e=0.48, Omega=0.0, i=30.0, w=0.0)
    orbit = pyasl.KeplerEllipse(a, per, e, Omega, i, w)
    t = np.linspace(0, 4, 200)
    pos = orbit.xyzPos(t)
    plt.plot(0, 0, 'bo', markersize=9, label="Earth")
    plt.plot(pos[::, 1], pos[::, 0], 'k-', label="Satellite Trajectory")
    plt.plot(pos[0, 1], pos[0, 0], 'r*', label="Periapsis")
    plt.legend(loc="upper right")
    plt.title('Orbital Simulation')
    plt.show()


    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D([0], [0], [0], 'bo',
            markersize=9, label="Earth")
    # Plot the satellite orbit
    ax.plot3D(pos[::, 1], pos[::, 0],
            pos[::, 2], 'k-',
            label="Satellite Trajectory")
    # Point of periapsis
    ax.plot3D([pos[0, 1]], [pos[0, 0]], [pos[0, 2]],
            'r*', label="Periapsis")
    plt.legend(loc="lower right")
    # Hide grid lines
    ax.grid(False)
    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.style.use('default')
    plt.title('Orbital Simulation')
    plt.show()



