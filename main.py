import math
from scipy.integrate import solve_ivp
import numpy as np
from density import *
from constants import *
from findpeaks import *
from us1976_analyses import *
import matplotlib.pyplot as plt
from ciraModel import *

# StudentModel or US1976 or CIRA or "Jacchia"

#method = "Jacchia"

#CD = 2.2
#m = 100
#A = math.pi * (1/2)**2
#R_a = 939 + R_earth #939
#R_p = 500 + R_earth #215
#Jacchia_temp = 600

# Read parameters from the file
with open('input.txt', 'r') as file:
    lines = file.readlines()

# Initialize variables
CD = None
m = None
A = None
R_a = None
R_p = None
method = None
Jacchia_temp = None

# Parse parameters
for line in lines:
    key_value = line.split(':')
    key = key_value[0].strip()
    
    if len(key_value) > 1:
        value = key_value[1].strip()
        if key == 'CD':
            # If the value is not there:
            if value == '':
                print("CD not specified. Please specify the parameters in the input.txt file.")
                exit()

            CD = float(value)
        elif key == 'm (Kg)':
            # If the value is not there:
            if value == '':
                print("m not specified. Please specify the parameters in the input.txt file.")
                exit()

            m = float(value)
        elif key == 'd (diameter m)':
            # If the value is not there:
            if value == '':
                print("d not specified. Please specify the parameters in the input.txt file.")
                exit()

            A = math.pi * (float(value)/2)**2
        elif key == 'Apoapsis Altitude (Km)':
            #If the value is not there:
            if value == '':
                print("R_a not specified. Please specify the parameters in the input.txt file.")
                exit()

            R_a = float(value) + R_earth
        elif key == 'Periapsis Altitude (Km)':
            #If the value is not there:
            if value == '':
                print("R_p not specified. Please specify the parameters in the input.txt file.")
                exit()

            R_p = float(value) + R_earth
        elif key == 'Method (StudentModel, US1976 , CIRA, Jacchia, MET)':
            # If the value is not there:
            if value == '':
                print("Method not specified. Please specify the parameters in the input.txt file.")
                exit()
            method = value
        elif key == 'Model Temperature (For Jacchia and MET) (Kelvin)':
            # If the value is not there:
            if value == '':
                print("Temperature not specified. Please specify the parameters in the input.txt file.")
                exit()
            Temp = float(value)

# Now, you have the parameters in the variables CD, m, A, R_a, R_p, method, and Jacchia_temp
#print("Parameters:")
#print(f"CD: {CD}")
#print(f"m (Kg): {m}")
#print(f"A (diameter m): {A}")
#print(f"Apoapsis Altitude (m): {R_a}")
#print(f"Periapsis Altitude (m): {R_p}")
#print(f"Method (StudentModel, US1976 , CIRA, Jacchia, MET): {method}")
#print(f"Model Temperature (For Jacchia and MET) (Kelvin): {Jacchia_temp}")

if method == "Jacchia":
    Jacchia_temp = Temp
elif method == "MET":
    MET_temp = Temp


# Eccentric anomaly
e = (R_a - R_p) / (R_a + R_p)

# Semi-major axis
a = (R_a + R_p) / 2

# Considering the satellite starts its journey in the apoapsis
R0 = R_a

# Function to detect the end of the integration
def event(t, y):
    altitude = np.linalg.norm(y[:2]) - R_earth

    if altitude <= 0:
        # if the altitude is less than or equal to 0, the integration stops
        print("Termination event triggered.")
        return 0
    return 1

def f(t, input):
    x, y, vx, vy = input
    r = (x, y)
    v = (vx, vy)
    r_unit = r / np.linalg.norm(r)
    r_norm = np.linalg.norm(r)
    r_perpend_unit = np.array([-r_unit[1], r_unit[0]])


    v_atm = omega_earth * r_norm

    vrel = v - float(v_atm) * r_perpend_unit
    
    vrel_abs = np.linalg.norm(vrel)
    vrel_unit = vrel / vrel_abs

    if method == "US1976":
        density = density_US1976(r_norm - R_earth)
    elif method == "StudentModel":
        density = student_model(r_norm - R_earth)
    elif method == "CIRA":
        density = density_CIRA(r_norm - R_earth)
    elif method == "Jacchia":
        density = density_jacchia(r_norm - R_earth, Jacchia_temp)
    elif method == "MET":
        density = density_MET(r_norm - R_earth, MET_temp)
    
    P = -CD* density * (1000 * vrel_abs)**2  * A  / (2 *m) * vrel_unit #m/s^2
    g = -mu / r_norm**2 * r_unit # Km/s^2
    a_fin = P / 1000 + g #Km/s^2

    return [vx, vy, a_fin[0], a_fin[1]]

def main():
    x0 = 0
    y0 = R0

    # vis-viva equation
    vx0 = math.sqrt(mu * (2 / (y0) - 1 / a)) #Km/s
    vy0 = 0

    t0 = 1e-10
    tf = 120 * days

    # Set error tolerances, initial step size, and termination event:
    options = {'rtol': 1e-7, 'atol': 1e-7}

    event.terminal = True

    sol = solve_ivp(f, [t0, tf], [x0, y0, vx0, vy0], method='DOP853', events=event)

    # Compute the altitudes in km in norm
    altitude = np.linalg.norm(sol.y[:2], axis=0) - R_earth

    time = sol.t * 10

    [i_max, i_min] = find_local_extrema(altitude, time)
    [max, min] = [altitude[i_max], altitude[i_min]]
    [tmax, tmin] = [time[i_max], time[i_min]]

    if method == "US1976":
        fig1 = plt.figure(1)

        # Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs
        degree = 10  # Degree of the polynomial curve to fit the data
        maxima_coefficients, maxima_fitted_values = fit_curve(tmax, max, degree)
        minima_coefficients, minima_fitted_values = fit_curve(tmin, min, degree)
        
        # Plotting maxima with fitted values
        plt.plot(tmax/86400, maxima_fitted_values, label='Maxima', linestyle='-', color='red')

        # Grid to the plot
        plt.grid()

        # Plotting minima with fitted values
        plt.plot(tmin/86400, minima_fitted_values, label='Minima', linestyle='-', color='green')
        plt.xlabel('Time (days)')
        plt.ylabel('Altitude (Km)')        

        # Display the plot
        #fig1.show()

        model_analyses(maxima_fitted_values, minima_fitted_values, tmax, tmin)

    if method == "CIRA":
        fig1 = plt.figure(1)

        # Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs
        degree = 10  # Degree of the polynomial curve to fit the data
        maxima_coefficients, maxima_fitted_values = fit_curve(tmax, max, degree)
        minima_coefficients, minima_fitted_values = fit_curve(tmin, min, degree)
        
        # Plotting maxima with fitted values
        plt.plot(tmax/86400, maxima_fitted_values, label='Maxima', linestyle='-', color='red')

        # Grid to the plot
        plt.grid()

        # Plotting minima with fitted values
        plt.plot(tmin/86400, minima_fitted_values, label='Minima', linestyle='-', color='green')
        plt.xlabel('Time (days)')
        plt.ylabel('Altitude (Km)')        

        # Display the plot
        #fig1.show()

        model_analyses(maxima_fitted_values, minima_fitted_values, tmax, tmin)

    elif method == "StudentModel":
        plt.plot(time/3600, altitude, label='Minima', linestyle='-', color='green')

        # Grid to the plot
        plt.grid()

        # Adding labels and legend

        plt.xlabel('Time (hours)')
        plt.ylabel('Altitude (Km)')

        plt.show()

    elif method == "Jacchia":
        fig1 = plt.figure(1)

        # Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs
        degree = 10  # Degree of the polynomial curve to fit the data
        maxima_coefficients, maxima_fitted_values = fit_curve(tmax, max, degree)
        minima_coefficients, minima_fitted_values = fit_curve(tmin, min, degree)
        
        # Plotting maxima with fitted values
        plt.plot(tmax/86400, maxima_fitted_values, label='Maxima', linestyle='-', color='red')

        # Grid to the plot
        plt.grid()

        # Plotting minima with fitted values
        plt.plot(tmin/86400, minima_fitted_values, label='Minima', linestyle='-', color='green')
        plt.xlabel('Time (days)')
        plt.ylabel('Altitude (Km)')        

        # Display the plot
        #fig1.show()

        model_analyses(maxima_fitted_values, minima_fitted_values, tmax, tmin)
    
    elif method == "MET":
        fig1 = plt.figure(1)

        # Fit the data aquiared to a polynomial curve in order to improve the visualization of the graphs
        degree = 10  # Degree of the polynomial curve to fit the data
        maxima_coefficients, maxima_fitted_values = fit_curve(tmax, max, degree)
        minima_coefficients, minima_fitted_values = fit_curve(tmin, min, degree)
        
        # Plotting maxima with fitted values
        plt.plot(tmax/86400, maxima_fitted_values, label='Maxima', linestyle='-', color='red')

        # Grid to the plot
        plt.grid()

        # Plotting minima with fitted values
        plt.plot(tmin/86400, minima_fitted_values, label='Minima', linestyle='-', color='green')
        plt.xlabel('Time (days)')
        plt.ylabel('Altitude (Km)')        

        # Display the plot
        #fig1.show()

        model_analyses(maxima_fitted_values, minima_fitted_values, tmax, tmin)

main()
