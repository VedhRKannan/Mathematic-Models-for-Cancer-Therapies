import numpy as np
import matplotlib.pyplot as plt

# Define Hill equation function


def hill(x, Kd, n):
    return (x**n)/(Kd + x**n)


# Define ligand concentration range
ligand = np.logspace(-3, 3, 1000)

# Define Kd and Hill coefficient
Kd = 1.0  # dissociation constant
n = 10  # Hill coefficient

# Compute receptor activation
receptor_activity = hill(ligand, Kd, n)

# Plot results
plt.semilogx(ligand, receptor_activity)
plt.xlabel('Ligand concentration (nM)')
plt.ylabel('Receptor activity')
plt.title('EGFR activation by EGF ligand binding')
plt.show()
