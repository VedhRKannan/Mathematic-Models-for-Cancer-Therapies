import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# import numpy as np

# Your data (substrate concentration and product concentration)
substrate_concentration = np.array([1, 10])
product_concentration = np.array([0.0564, 0.1535])

# Incubation time
# in minutes (example value, replace with your actual time)
incubation_time = 20

# Calculate reaction rate (product concentration / incubation time)
reaction_rate = product_concentration / incubation_time


def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


# Fit the curve using SciPy's curve_fit function
params, cov = curve_fit(
    michaelis_menten, substrate_concentration, reaction_rate)

# Extract fitted parameters Vmax and Km
Vmax, Km = params


# Generate a range of substrate concentrations for the fitted curve
substrate_range = np.linspace(0, max(substrate_concentration), 100)

# Calculate the fitted reaction rates for the substrate range using the fitted parameters
fitted_reaction_rate = michaelis_menten(substrate_range, Vmax, Km)
print(Km)
print(Vmax)
# plt.plot(substrate_concentration, reaction_rate,
#  'bo', label='Experimental data')
plt.plot(substrate_range, fitted_reaction_rate,
         'r-', label='Fitted curve')
plt.grid()
plt.xlabel('Substrate concentration (mM)')
plt.ylabel('Reaction rate (mM/min)')
plt.title('Michaelis-Menten Model')
plt.axvline(x=2.3654, label='Km')
plt.legend()
plt.show()
