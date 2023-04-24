import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the hill model function


def hill_model(x, A, n, Kd):
    return A * (x ** n) / (Kd ** n + x ** n)


# Define the experimental data
EGF_concentration = np.array([0,])  # EGF concentration in nM
# Fluorescence intensity in arbitrary units
fluorescence_intensity = np.array([0,])

# Define initial guess for fitting parameters
initial_guess = [100, 1, 1]

# Fit the experimental data to the hill model
popt, pcov = curve_fit(hill_model, EGF_concentration,
                       fluorescence_intensity, p0=initial_guess)

# Extract the fitting parameters and calculate the hill coefficient
A_fit, n_fit, Kd_fit = popt
n_hill = n_fit

# Print the hill coefficient
print("The hill coefficient (n) is:", n_hill)

# Plot the experimental data and the fitted curve
plt.plot(EGF_concentration, fluorescence_intensity,
         'ro', label='Experimental data')
plt.plot(EGF_concentration, hill_model(
    EGF_concentration, A_fit, n_fit, Kd_fit), 'b-', label='Fit')
plt.legend()
plt.xlabel('EGF concentration (nM)')
plt.ylabel('Fluorescence intensity (AU)')
plt.show()
