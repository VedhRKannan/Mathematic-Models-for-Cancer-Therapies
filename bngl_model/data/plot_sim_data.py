import numpy as np
import matplotlib.pyplot as plt

cdat_file = "bngl_model/data/egfr_growth_normal.cdat"
# cdat_file = "bngl_model/data/egfr_growth_overexpressed.cdat"


with open(cdat_file, "r") as f:
    lines = f.readlines()

data = []
for line in lines:
    if line.startswith("#"):
        continue
    values = line.strip().split()
    data.append([float(x) for x in values])

data = np.array(data)

time = data[:, 0]
EGFR_monomer = data[:, 4]
EGF = data[:, 1]
EGFR_dimer = data[:, 7]

# print(EGFR_dimer)
print(time)
# print(EGF)
# print(EGFR_monomer)
fig, ax = plt.subplots()
ax.plot(time, EGFR_monomer, label="EGFR monomer")
ax.plot(time, EGF, label="EGF")
ax.plot(time, EGFR_dimer, label="EGFR dimer")
ax.set_yscale('log', base=10)
ax.set_xscale('log', base=10)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.title("EGFR Dimer Growth")
plt.legend()

# plt.savefig("EGFR_dimer_growth_python.png")
plt.show()
