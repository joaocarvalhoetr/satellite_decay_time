from constants import *
import matplotlib.pyplot as plt
from findpeaks import *

def model_analyses(maxima_fitted_values, minima_fitted_values, tmax, tmin):

    l_max = len(maxima_fitted_values)
    l_min = len(minima_fitted_values)

    e_values = []
    a_values = []
    time_list = []

    # Compare the lenght of max and min and remove the excess of the longest one
    if l_max > l_min:
        maxima_fitted_values = maxima_fitted_values[:l_min]
        tmax = tmax[:l_min]
    elif l_min > l_max:
        minima_fitted_values = minima_fitted_values[:l_max]
        tmin = tmin[:l_max]

    for i in range(len(maxima_fitted_values)):
        R_a = maxima_fitted_values[i] + R_earth
        R_p = minima_fitted_values[i] + R_earth

        time_list.append(tmax[i]/86400)

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

    


