from constants import *
import matplotlib.pyplot as plt
from findpeaks import *

def model_analyses(i_max, i_min, altitude, time):

    l_max = len(i_max)
    l_min = len(i_min)

    e_values = []
    a_values = []
    time_list = []

    # Compare the lenght of max and min and remove the excess of the longest one
    if l_max > l_min:
        i_max = i_max[:l_min]
    elif l_min > l_max:
        i_min = i_min[:l_max]

    for i in range(len(i_min)):
        R_a = altitude[int(i_max[i])] + R_earth
        R_p = altitude[int(i_min[i])] + R_earth

        time_list.append(time[int(i_max[i])]/86400)

        e = (R_a - R_p) / (R_a + R_p)
        e_values.append(e)

        # Semi-major axis
        a = (R_a + R_p) / 2
        a_values.append(a)

    # Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs
    degree = 10  # Degree of the polynomial curve to fit the data
    maxima_coefficients, e_fitted_values = fit_curve(time_list, e_values, degree)
    minima_coefficients, a_fitted_values = fit_curve(time_list, a_values, degree)


    
    # Plot the eccentricity for the times and the semi-major axis
    fig2 = plt.figure(2)
    plt.plot(time_list, e_fitted_values, label='Minima', linestyle='-', color='green')

    # Grid to the plot
    plt.grid()

    # Adding labels and legend

    plt.xlabel('Time (days)')
    plt.ylabel('Excentricity')

    #fig2.show()

    # Plot the eccentricity for the times and the semi-major axis
    fig3 = plt.figure(3)
    plt.plot(time_list, a_fitted_values, label='Minima', linestyle='-', color='red')

    # Grid to the plot
    plt.grid()

    # Adding labels and legend

    plt.xlabel('Time (days)')
    plt.ylabel('Semi-major Axis (Km)')

    plt.show()

    


