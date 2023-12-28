import matplotlib.pyplot as plt

# Your data
data = [253.39883203, 220.04789931, 218.03417112, 247.4292954, 304.07797653,
        383.01568924, 468.67966353, 552.87868285, 634.21152941, 701.26519296,
        748.72733358, 777.86206988, 782.2996517, 756.81017073, 709.17713866,
        642.38798563, 553.79297715, 451.19663837, 346.30993769, 248.22325343,
        166.21300097, 108.6204458]

# Plotting
plt.plot(data, marker='o', linestyle='-')

# Adding labels and title
plt.xlabel('X-axis')
plt.ylabel('Values')
plt.title('Plot of the given data')

# Display the plot
plt.show()
