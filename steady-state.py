import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-10, 10, 100)

y = np.zeros_like(x)

# plt.plot(x, y)

# plt.xlabel('x')
# plt.ylabel('d[ES]/dt')
# plt.title('Graph of d[ES]/dt = 0')

# plt.show()


Vmax = 1.0
Km = 0.5

# Define the substrate concentration range
S = np.linspace(0, 2, 100)

# Calculate the [ES] concentration as a function of [S]
ES = (S * Vmax) / (Km + S)

# Define the horizontal line for d[ES]/dt = 0
dES_dt = np.zeros_like(S)

# Plot the Michaelis-Menten curve and the horizontal line
plt.plot(S, ES, label='[ES]')
plt.plot(S, dES_dt, label='d[ES]/dt = 0')


plt.xlabel('[S]')
plt.ylabel('[ES]')
plt.title('Michaelis-Menten curve with d[ES]/dt = 0')
plt.legend()
plt.show()

# Add a legend
