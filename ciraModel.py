import numpy as np
import matplotlib.pyplot as plt
   

def interpolate_values(h_values, r_values, given_h):
        inter_r = np.interp(given_h, h_values, r_values)
        return inter_r

def density_CIRA(z):
    # Geometric altitudes (km):
    h = np.array([0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,
                  380,400,420,440,460,480,500,520,540,560,580
                  ,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900])

    # Corresponding densities (kg/m^3) from USSA76:
    r = np.array([1.16E+00,9.37E-02,4.02E-03,3.26E-04, 1.83E-05,5.73E-07,2.03E-08,3.44E-09,1.20E-09,5.46E-10,
                  2.84E-10,1.61E-10,	9.60E-11,
                  5.97E-11,3.83E-11,2.52E-11,1.69E-11,1.16E-11,7.99E-12,5.60E-12,3.96E-12,2.83E-12,2.03E-12,
                  1.47E-12,1.07E-12,7.85E-13,5.78E-13,4.29E-13,3.19E-13,2.39E-13,1.80E-13,1.36E-13,1.04E-13,
                  7.98E-14,6.16E-14,4.80E-14,3.76E-14,2.98E-14,2.38E-14,1.92E-14,1.57E-14,1.29E-14,1.07E-14,9.03E-15,7.67E-15,6.59E-15])

  

    if z <=180 :
        density=interpolate_values(h,r,z)
        print(density)
        return density
    
    if z>180 or z<=340 :
        # Select the range of data for fitting
        fit_range_start = 180
        fit_range_end = 340
        mask = (h >= fit_range_start) & (h <= fit_range_end)
        fit_h = h[mask]
        fit_r = r[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 2

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)
        
        density = poly_function(z)
        
        return density

    if z>340 or z<=480 :
        # Select the range of data for fitting
        fit_range_start = 340
        fit_range_end = 480
        mask = (h >= fit_range_start) & (h <= fit_range_end)
        fit_h = h[mask]
        fit_r = r[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 2

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)

        density = poly_function(z)
        return density
    
    if z>480 or z<=640 :
        # Select the range of data for fitting
        fit_range_start = 480
        fit_range_end = 640
        mask = (h >= fit_range_start) & (h <= fit_range_end)
        fit_h = h[mask]
        fit_r = r[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 2

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)

        density = poly_function(z)
        return density

    if z>640 or z<=840 :
        # Select the range of data for fitting
        fit_range_start = 640
        fit_range_end = 840
        mask = (h >= fit_range_start) & (h <= fit_range_end)
        fit_h = h[mask]
        fit_r = r[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 2

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)

        density = poly_function(z)
        return density

    if z>840 or z<=900 :
        # Select the range of data for fitting
        fit_range_start = 840
        fit_range_end = 900
        mask = (h >= fit_range_start) & (h <= fit_range_end)
        fit_h = h[mask]
        fit_r = r[mask] 

        # Choose the degree of the polynomial (2 in this case)
        degree = 2

        # Fit a polynomial to the data
        coefficients = np.polyfit(fit_h, fit_r, degree)

        # Create a polynomial function based on the coefficients
        poly_function = np.poly1d(coefficients)

        density = poly_function(z)
        return density

    if z > 900:
         return 6.59E-15

