import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("EGFR_dimer_growth.txt", delimiter=",", skiprows=1)

time = data[:, 0]
EGFR_monomer = data[:, 1]
EGF = data[:, 2]
EGFR_dimer = data[:, 3]

plt.plot(time, EGFR_monomer, label="EGFR monomer")
plt.plot(time, EGF, label="EGF")
plt.plot(time, EGFR_dimer, label="EGFR dimer")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.title("EGFR Dimer Growth")
plt.legend()

plt.savefig("EGFR_dimer_growth_python.pdf")
plt.show()
