import math
from constants import *
import pandas as pd
from io import StringIO
import numpy as np
import os

def line_find(height):
    file = pd.read_csv('atmosphericlayers_output.csv', sep=',', header=0)
    for i in range(20):
        i = i+1
        #find the layer for the given height
        current_line = file.iloc[i]
        Z = current_line['Geometric Altitude $Z$ (km)']
        if Z > height:
            i = i-1
            return file.iloc[i]


# Compute the densities for the atmospheric layers
file = pd.read_csv('atmosphericlayers.csv', sep=',', header=0)
file.to_csv('atmosphericlayers_output.csv', index=False)
for i in range(20):
    i = i+1

    file = pd.read_csv('atmosphericlayers_output.csv', sep=',', header=0)
    current_line = file.iloc[i]

    # In the current layer
    Z = current_line['Geometric Altitude $Z$ (km)']
    T_M = current_line['Molecular Temperature $T_M$ (K)']

    #Read the density_0 and Z_bef from the line before
    line_before = file.iloc[i-1]
    density_0 = line_before['density_0']
    
    Z_bef = line_before['Geometric Altitude $Z$ (km)']

    density = density_0*math.exp(-g_0*(Z-Z_bef)/(R_idealgas*T_M))

    # insert density in the current_line in the column density_0
    file.at[i, 'density_0'] = density

    #save the new file in the csv
    file.to_csv('atmosphericlayers_output.csv', index=False)



# Function that calculates a specific temperature for a given height. It uses the data from the csv file
# in order to with the lapse rate compute the temperature wanted.
  
def temperature_calculate(heigh):
    
    # Locate the layer the height given is in
    current_line = line_find(heigh)

    if current_line is None:
        current_line = line_find(2600)
    Z = current_line['Geometric Altitude $Z$ (km)']

    #Compute the temperature for the given height
    temperature = current_line['Molecular Temperature $T_M$ (K)'] + current_line['Lapse Rate $L_h$ (K/km)']*(heigh-Z)

    #Return the temperature
    return temperature

# Function that calculates the scale factor for a given height. It uses the function temperature_calculate as auxiliary.

def scale_factor_calculate(height):
    
    #Compute the temperature for the given height
    temperature = temperature_calculate(height)
    #Compute the scale factor for the given height
    scale_factor = (R_idealgas*temperature)/(g_0)
    #Return the scale factor
    return scale_factor

# Function that calculates the density for a given height. It uses the function scale_factor_calculate as auxiliary.

def density_calculate(heigh):
    #Compute the scale factor for the given height
    scale_factor = scale_factor_calculate(heigh)

    # Locate the layer the height given is in
    current_line = line_find(heigh)

    if current_line is None:
        current_line = line_find(2600)
    Z = current_line['Geometric Altitude $Z$ (km)']

    #Read the density_0 from the line
    density_0 = current_line['density_0']

    #Compute the density for the given height
    density = density_0*math.exp(-heigh/scale_factor)

    #Return the density
    return density

print(density_calculate(254))