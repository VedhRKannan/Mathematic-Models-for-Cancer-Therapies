from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# Define the Michaelis-Menten equation


def michaelis_menten(x, vmax, km):
    return (vmax * x) / (km + x)

# Define a function to calculate the effective enzyme concentration based on pH and temperature


def effective_enzyme_concentration(enzyme_concentration, ph, temperature, denaturation_constant):
    effective_concentration = enzyme_concentration * \
        10**((ph - 7.0) * denaturation_constant)
    effective_concentration *= np.exp(-denaturation_constant *
                                      (temperature - 37.0))
    return effective_concentration

# Define the reaction function


def reaction(t, y, vmax, km, enzyme_concentration, ph, temperature, denaturation_constant):
    substrate_concentration = y[0]
    enzyme_effective_concentration = effective_enzyme_concentration(
        enzyme_concentration, ph, temperature, denaturation_constant)
    velocity = michaelis_menten(
        substrate_concentration, vmax, km) * enzyme_effective_concentration
    dydt = -velocity, velocity
    return dydt


# Define the initial conditions
y0 = [1.0, 0.0]

# Define the parameters
vmax = 10.0
km = 0.1
enzyme_concentration = 10
ph = 8.6
temperature = 37.0
denaturation_constant = 0.01

# Define the time span and step size
t_span = [0.0, 10.0]
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Solve the differential equation
solution = solve_ivp(reaction, t_span, y0, args=(
    vmax, km, enzyme_concentration, ph, temperature, denaturation_constant), t_eval=t_eval)

# Plot the results
plt.plot(solution.t, solution.y[0], label='substrate')
plt.plot(solution.t, solution.y[1], label='product')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.show()
