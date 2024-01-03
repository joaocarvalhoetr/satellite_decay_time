import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Given data
time = np.array([0, 10, 20, 30, 500])
altitude = np.array([500, 400, 300, 200, 0])  # Replace ... with actual values
velocity_x = np.array([70, 100, 200, 300, 400])  # Replace ... with actual values

# Orbital parameters
earth_radius = 6371  # Earth radius in kilometers

# Calculate orbital position
theta = np.radians(time)  # Convert time to radians for equatorial plane

# Calculate x and y coordinates in the equatorial plane
x = (earth_radius + altitude) * np.cos(theta)
y = (earth_radius + altitude) * np.sin(theta)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the orbit
ax.plot(x, y, altitude, label='Orbit')

# Set labels
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Altitude (km)')
ax.set_title('3D Earth Orbit')

# Show the plot
plt.show()
