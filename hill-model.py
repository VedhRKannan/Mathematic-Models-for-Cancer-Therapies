import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def hill_equation(x, bmax, kd, n):
    return bmax * (x ** n) / ((kd ** n) + (x ** n))


egf_conc = np.array([0.1, 0.3, 1, 3, 10])   # EGF concentrations in nM
# measured binding values in fmol/mg protein
binding = np.array([2.1, 5.8, 16.3, 46.7, 88.6])

max_binding = np.max(binding)   # maximum binding capacity
binding_norm = binding / max_binding   # normalize to maximum binding
binding_norm = binding_norm * 1000   # convert to fmol/μg protein

p0 = [100, 1e-9, 2]   # initial guesses for the parameters
popt, pcov = curve_fit(hill_equation, egf_conc, binding_norm, p0)
bmax_fit, kd_fit, n_fit = popt

plt.plot(egf_conc, binding_norm, 'ro', label='Experimental Data')
plt.plot(egf_conc, hill_equation(egf_conc, bmax_fit,
         kd_fit, n_fit), 'b-', label='Hill Equation Fit')
plt.legend()
plt.xlabel('EGF Concentration (nM)')
plt.ylabel('Normalized Binding (fmol/μg protein)')
plt.show()
