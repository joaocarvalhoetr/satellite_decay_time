import math

# Constants with units
mu = 398600.4418  # Earth's gravitational parameter in km^3/s^2
R_earth = 6378.137  # Earth's radius in km
hours = 3600  # seconds in an hour
days = 24 * hours  # seconds in a day
R_idealgas = 8.31446261815324  # Ideal gas constant in J/(mol K)
R_air = 287.05  # gas constant for air
g_0 = 9.80665  # Acceleration due to gravity at Earth's surface in m/s^2
rho_0 = 1.225  # Standard sea-level atmospheric density in kg/m^3
deg = math.pi / 180  # Conversion factor from degrees to radians
omega_earth = 2 * math.pi / (24 * 3600)  # Earth's angular velocity in rad/s
wE = 7.2921159e-5  # Earth's angular velocity vector in rad/s
