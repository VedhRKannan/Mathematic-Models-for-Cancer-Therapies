import numpy as np
import matplotlib.pyplot as plt

cdat_file = "bngl_model/histidine/2023_11_13__14_53_52/histidine.cdat"
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
Tyrosine = data[:, 2]
Histidine = data[:, 1]
Phosphate = data[:, 3]
Phosphohistidine = data[:, 4]
Phosphotyrosine = data[:, 5]

# print(EGFR_dimer)
print(Histidine)
# print(EGF)
# print(EGFR_monomer)
fig, ax = plt.subplots()
ax.plot(time, Tyrosine, label="Tyrosine")
ax.plot(time, Histidine, label="Histidine")
ax.plot(time, Phosphate, label="Phosphate")
ax.plot(time, Phosphohistidine, label="Phosphohistidine")
ax.plot(time, Phosphotyrosine, label="Phosphotyrosine")
ax.set_yscale('log', base=10)
ax.set_xscale('log', base=10)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.title("EGFR Dimer Growth")
plt.legend()

# plt.savefig("EGFR_dimer_growth_python.png")
plt.show()
