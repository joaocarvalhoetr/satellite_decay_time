import math
from scipy.integrate import solve_ivp
import numpy as np
from density import *
from constants import *
from findpeaks import *
from analyses import *
import matplotlib.pyplot as plt
from ciraModel import *

# Header Comments
#
# Project: Satellite Reentry
#
# Description: This program computes the time of a satellite into the Earth's atmosphere. It solves an ODE system that
# describes the satellite's motion. The program also computes the altitude of the satellite as a function of time and
# plots the results. For each method, the program computes the specific density for each time step.
# The ODE system is solved using the DOP853 which is a Runge-Kutta method of order 8. There's a TXT file that contains
# the parameters of the satellite and the method to be used. The program reads the parameters from the file and computes
# the results.
# There's 5 methods available: StudentModel, US1976, CIRA, Jacchia, and MET. The StudentModel is a model developed by
# the students that considers g constant and uses Scale Factors.
#
# We took inspiration from the book: Howard D. Curtis - Orbital Mechanics For Engineering Students
#
# Authors: 
# António José Domingos Reis (102473), Ricardo Gandra de Sousa (102498), Tiago André da Silva Ruge (102551), 
# João P. Veloso Onofre de Carvalho (102686), Eduardo De Almeida Helena (102793), Tomás Gomes Coelho (102805),
# Gonçalo José Reis Bessa da Silva (102995), Fernando Meneses Vicente (103048), Álvaro G. Silva Vilela Caridade (103526)
# João Nuno Rodrigues Alegrete (103676), João Bessa e Silva Machado Vilaça (103966)
#
# Current Version: 1.0 (Jan 5, 2024)
#
#
# Version History:
# Version 1.0 (Jan 5, 2024)
# - Updated temperature values in models.
# - Renamed analysis-related components.
# - Merged changes from the main branch.
# - Added MET (Mission Elapsed Time) to the main branch.
# - Removed obsolete code and files.
# - Work related to MET density.

# Version 0.9 (Jan 4, 2024)
# - Updated density calculations - CIRA, Nasa Model, Student Developed.

# Version 0.8 (Jan 2, 2024)
# - Updated constants in the code.
# - Changed the name of the acceleration function.

# Version 0.7 (Jan 1, 2024)
# - Organized files, moving them to the TEMP directory.

# Version 0.6 (Dec 30, 2023)
# - Deleted altitude_vs_time.gif.
# - Deleted 2D and 3D orbit graphics.

# Version 0.5 (Dec 29, 2023)
# - Added 2D and 3D orbit graphics as the satellite falls.
# - Rotated the Earth in the graphics.

# Version 0.4 (Dec 28, 2023)
# - Created function to mask data.
# - Updated the main script.

# Version 0.3 (Dec 27, 2023)
# - Main modifications.
# - Updated simplified model.
# - Updated the simplified test script.

# Version 0.2 (Dec 26, 2023)
# - Introduced a simplified version of the program.

# Version 0.1 (Dec 23, 2023)
# - Added 3D animation functionality.

# Version 0.0 (Dec 17, 2023)
# - ODE (Ordinary Differential Equation) solver is working.
# - Work on the ODE system.

# Initial Commit (Dec 16, 2023)
# - Initial working version.
# - Creation of US1976 Model.
# - Added functionality to read from text files.
#
# End of Header Comments

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
                input("Press Enter to close the terminal...")
                exit()

            CD = float(value)
        elif key == 'm (Kg)':
            # If the value is not there:
            if value == '':
                print("m not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()

            m = float(value)
        elif key == 'd (diameter m)':
            # If the value is not there:
            if value == '':
                print("d not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()

            A = math.pi * (float(value)/2)**2
        elif key == 'Apoapsis Altitude (Km)':
            #If the value is not there:
            if value == '':
                print("R_a not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()

            R_a = float(value) + R_earth
        elif key == 'Periapsis Altitude (Km)':
            #If the value is not there:
            if value == '':
                print("R_p not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()

            R_p = float(value) + R_earth
        elif key == 'Method (StudentModel, US1976 , CIRA, Jacchia, MET)':
            # If the value is not there:
            if value == '':
                print("Method not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()
            method = value
        elif key == 'Model Temperature (For Jacchia and MET) (Kelvin)':
            # If the value is not there:
            if value == '':
                print("Temperature not specified. Please specify the parameters in the input.txt file.")
                input("Press Enter to close the terminal...")
                exit()
            # Check if the temperature is available Jaccia or MET. Jaccia available is 600 1500 2400 and MET 600 1400 and 2200
            
            if method == "Jacchia":
                if value == '600' or value == '1500' or value == '2400':
                    pass
                else:
                    print("Temperature not available for Jaccia. Please specify the parameters in the input.txt file.")
                    input("Press Enter to close the terminal...")
                    exit()

            elif method == "MET":
                if value == '600' or value == '1400' or value == '2200':
                    pass
                else:
                    print("Temperature not available for MET. Please specify the parameters in the input.txt file.")
                    input("Press Enter to close the terminal...")
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
