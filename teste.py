import numpy as np
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import PyAstronomy.pyasl as pyasl
from PyAstronomy import pyasl



orbit = pyasl.KeplerEllipse(a=1, per=1, e=0.5, Omega=0.0, i=30.0, w=0.0)
#orbit = pyasl.KeplerEllipse(a, per, e, Omega, i, w)
t = np.linspace(0, 40, 200)
pos = orbit.xyzPos(t)




fig = plt.figure()
ax = p3.Axes3D(fig)
red_dot, = ax.plot(pos[::, 1], pos[::, 0], pos[::, 2], 'ro',label="Satellite")


def animate(i, pos, red_dot):
    #red_dot.set_data([pos[i][1], pos[i][0]])
    #red_dot.set_3d_properties(pos[i][2])
    red_dot.set_data(pos[:i, 1], pos[:i, 0])
    red_dot.set_3d_properties(pos[:i, 2])
 
    return red_dot,

# create animation using the animate() function
ani = animation.FuncAnimation(fig, animate, frames=len(pos), fargs=(pos, red_dot), interval=100,blit=False)

#ax.set_xlim([np.min(pos[:, 1]), np.max(pos[:, 1])])
#ax.set_ylim([np.min(pos[:, 0]), np.max(pos[:, 0])])
#ax.set_zlim([np.min(pos[:, 2]), np.max(pos[:, 2])])



ax.plot([0], [0], [0], 'bo', markersize=9, label="Earth")
ax.plot(pos[::, 1], pos[::, 0], pos[::, 2], 'k-',label="Satellite Trajectory")
ax.legend()

# Hide grid lines
ax.grid(False)
# Hide axes ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
plt.style.use('default')
plt.legend()
plt.show()