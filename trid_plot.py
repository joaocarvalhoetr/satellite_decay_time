import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def tridimensionalPlot(x,y) :
    # Generate 2D data
    #x = np.linspace(-5, 5, 100)
    #y = np.sin(x)

    # Create a 3D plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Extrude the 2D plot into 3D
    z = np.zeros_like(x)  # Set the z-values to a constant (e.g., 0)
    ax.plot(x, y, z, color='b', label='2D Plot in 3D')

    # Customize the plot
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    ax.legend()

    # Show the plot
    plt.show()
